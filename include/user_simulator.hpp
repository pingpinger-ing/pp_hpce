
#ifndef simulator_hpp
#define simulator_hpp

#include <cstdint>
#include <vector>
#include <deque>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include <stdarg.h>

#include "util.hpp"

#include <tbb/parallel_for.h>

template<class TGraph>
class Simulator
{
public:
    typedef typename TGraph::graph_type graph_type;
    typedef typename TGraph::properties_type properties_type;
    typedef typename TGraph::device_type device_type;
    typedef typename TGraph::message_type message_type;
    typedef typename TGraph::channel_type channel_type;
    typedef typename TGraph::SupervisorDevice SupervisorDevice;
private:    
    struct node;
    struct edge;
    struct output;

    struct node
    {
        properties_type properties;
        device_type state;
        std::vector<edge*> incoming;
        std::vector<edge*> outgoing;
        bool output;
    };
    
    struct edge
    {
        node *dst;
        node *src;
        unsigned delay;         // How long it takes a message to get through
        
        channel_type channel;
        
        unsigned messageStatus; // 0->empty, 1->ready, 2->inflight
        message_type messageData;
    };
    
    struct output
    {
        const properties_type *source;  // Where the output came from
        message_type output;            // Message associated with the output
        unsigned sendStep;              // Which step was it send in?
    };
    
    struct stats
    {
        uint32_t stepIndex;
        
        uint32_t nodeIdleSteps;
        uint32_t nodeBlockedSteps;
        uint32_t nodeSendSteps;
        
        uint32_t edgeIdleSteps;
        uint32_t edgeTransitSteps;
        uint32_t edgeDeliverSteps;
    };
    
    int m_logLevel;
    
    void log(int level, const char *msg, ...)
    {
        if(level <= m_logLevel){
            char localBuffer[256];
            char *globalBuffer=0;
            char *buffer=localBuffer;
            
            va_list va;
            va_start(va, msg);
            int n=vsnprintf(buffer, sizeof(localBuffer), msg, va);
            va_end(va);
            
            if(n<=0){
                throw std::runtime_error("log failure.");
            }
            
            if(n >= (int)sizeof(localBuffer)){
                globalBuffer=new char[n+1];
                buffer=globalBuffer;
                va_list va;
                va_start(va, msg);
                vsnprintf(buffer, n+1, msg, va);
                va_end(va);
            }
            
            
            
            fprintf(stderr, "[Sim], %u, %.3f, %s\n", level, puzzler::now()*1e-9, buffer);
            
            if(globalBuffer){
                delete []globalBuffer;
            }
        }
    }
    
    uint32_t m_step;
    graph_type m_graph;
    std::vector<node> m_nodes;
    std::vector<edge> m_edges;
    std::deque<output> m_outputs;
    SupervisorDevice m_supervisor;
    
    std::ostream &m_statsDst;
    stats m_stats;
    
    // Give a single node (i.e. a device) the chance to
    // send a message.
    // \retval Return true if the device is blocked or sends. False if it is idle.
    
     bool update_node(unsigned index, node *n)
    {
        bool act = false;
        for (unsigned i = 0; i != n->incoming.size(); i++) {
            edge *e = n->incoming[i];
            switch (e->messageStatus) {
            case 0:
                continue;
            case 1:                // Deliver the message to the device                
                TGraph::on_recv(
                    &m_graph,
                    &(e->channel),
                    &(e->messageData),
                    &(n->properties),
                    &(n->state)
                );
            default:
                e->messageStatus--;
                act = true;
                continue;
            }
        }
        return act;
    }
    
    
    uint32_t stats_node(node *n)  // unsigned int 
    {
        if(!TGraph::ready_to_send(&m_graph, &(n->properties), &(n->state)) ){
            log(4, "  node %u : idle", index);
            m_stats.nodeIdleSteps++;
            return false; // Device doesn't want to send
        }
        
    for(unsigned i=0; i < n->outgoing.size(); i++){
            if( n->outgoing[i]->messageStatus>0 ){
                m_stats.nodeBlockedSteps++;
                return true; // One of the outputs is full, so we are blocked
            }
        }
     message_type message;
    
    // Get the device to send the message
        n->output = TGraph::on_send(
            &m_graph,
            &message,
            &(n->properties),
            &(n->state)
        );
        
        for(unsigned i=0; i < n->outgoing.size(); i++){
            assert( 0 == n->outgoing[i]->messageStatus );
            n->outgoing[i]->messageData = message; // Copy message into channel
            n->outgoing[i]->messageStatus = 1 + n->outgoing[i]->delay; // How long until it is ready?
        }

        m_stats.nodeSendSteps++;
        return true;
    }
    
    bool stats_edge(const edge *e)
    {
        switch (e->messageStatus) {
        case 0:
            m_stats.edgeIdleSteps++;
            return false;
        case 1:
            m_stats.edgeDeliverSteps++;
            return true;
        default:
            m_stats.edgeTransitSteps++;
            return true;
        }
    }
    
     void stats_nodes(node *n, unsigned cnt, unsigned *idle, unsigned *blocked, unsigned *send)
    {
            std::vector<unsigned> p_idle(cnt), p_blocked(cnt), p_send(cnt);
            tbb::parallel_for(0u, cnt, [=, &p_idle, &p_blocked, &p_send](unsigned i) {
                stats_nodes(n, cnt, &p_idle[i], &p_blocked[i], &p_send[i]);
            });
 
    }
    
     bool stats_nodes()
    {
        unsigned idle, blocked, send;
        stats_nodes(m_nodes.data(), m_nodes.size(), &idle, &blocked, &send);
        m_stats.nodeIdleSteps += idle;
        m_stats.nodeBlockedSteps += blocked;
        m_stats.nodeSendSteps += send;
        return blocked || send;
    }
    
      
    bool step_all()
    {
          log(2, "stepping edges");
        bool active=false;
        // Edge statistics
        for (const edge &e: m_edges)
            active |= stats_edge(&e);
        
        
        tbb::parallel_for(tbb::blocked_range<unsigned>(0, m_nodes.size(), 512), [&](const tbb::blocked_range<unsigned>& range) {
            unsigned s = range.begin(), e = range.end();
            for (unsigned i = s; i != e; i++)
                update_node(i, &m_nodes[i]);
        }, tbb::simple_partitioner());
        
      log(2, "stepping nodes");
        // Node statistics
        active |= stats_nodes();  
        
       for (node &n: m_nodes)
            if (n.output) {  
                m_supervisor.onDeviceOutput(&(n.properties), &n.outgoing[0]->messageData);
                n.output = false;
            }
        
        return active;
    }
    
    
    
    void reset()
    {
        log(2, "resetting nodes");
        m_step=0;
        for(unsigned i=0; i<m_nodes.size(); i++){
            TGraph::on_init(&m_graph, &m_nodes[i].properties, &m_nodes[i].state);
        }
        log(2, "resetting edges");
        for(unsigned i=0; i<m_edges.size(); i++){
            m_edges[i].messageStatus = 0;
        }
    }
    
public:
    Simulator(
        int logLevel,
        std::ostream &stats,
        FILE *destFile,
        const graph_type &graph,
        unsigned numDevices, 
        unsigned numChannels
    )
        : m_logLevel(logLevel)
        , m_step(0)
        , m_graph(graph)
        , m_supervisor(&m_graph, destFile)
        , m_statsDst(stats)
    {
        m_nodes.reserve(numDevices);
        m_edges.reserve(numChannels);
    }
    
    
    unsigned addDevice(
        const properties_type &device
    ){
        unsigned index=m_nodes.size();
        node n;
        n.properties=device;
        m_nodes.push_back(n);
        
        m_supervisor.onAttachNode(&m_nodes[index].properties);
        
        return index;
    }
    
    void addChannel(
        unsigned srcIndex,
        unsigned dstIndex,
        unsigned delay,
        const channel_type &channel
    ){
        unsigned edgeIndex = m_edges.size();
        edge e;
        e.src = &m_nodes.at(srcIndex);
        e.dst = &m_nodes.at(dstIndex);
        e.delay = delay;
        e.channel = channel;
        e.messageStatus=0;
        m_edges.push_back(e);
        
        m_nodes.at(srcIndex).outgoing.push_back( &m_edges[edgeIndex] );
        m_nodes.at(dstIndex).incoming.push_back( &m_edges[edgeIndex] );
    }
    
    
    

    void run()
    {
        log(1, "begin run");
        
        bool active=true;
        
        reset();
        
        while(active){
            log(1, "step %u", m_step);
            
            m_stats={m_step, 0,0,0, 0,0,0};

            // Run all the nodes
            active = step_all();
            
            // Flush any outputs from the queue to the supervisor
            while(!m_outputs.empty()){
                const output &o = m_outputs.front();
                m_supervisor.onDeviceOutput(o.source, &o.output);
                m_outputs.pop_front();
            }
            
            // Send statistics out
            m_statsDst<<m_stats.stepIndex<<", "<<m_stats.nodeIdleSteps<<", "<<m_stats.nodeBlockedSteps<<", "<<m_stats.nodeSendSteps;
            m_statsDst<<", "<<m_stats.edgeIdleSteps<<", "<<m_stats.edgeTransitSteps<<", "<<m_stats.edgeDeliverSteps<<"\n";

            m_step++;            
         }
    }
};

#endif


