
#include "bdsg/snarl_distance_index.hpp"
#include "bdsg/ch.hpp" // BINDER_IGNORE because ch.hpp is C++20 + Boost and is not bound
#include <jansson.h>
#include <arpa/inet.h>

using namespace std;
using namespace handlegraph;
namespace bdsg {

const vector<std::string> SnarlDistanceIndex::record_t_as_string = {
    "ROOT", "NODE", "DISTANCED_NODE", 
    "TRIVIAL_SNARL", "DISTANCED_TRIVIAL_SNARL",
    "SIMPLE_SNARL", "DISTANCED_SIMPLE_SNARL",
    "REGULAR_SNARL", "DISTANCED_REGULAR_SNARL", "OVERSIZED_REGULAR_SNARL",
    "SNARL", "DISTANCED_SNARL", "OVERSIZED_SNARL", 
    "ROOT_SNARL", "DISTANCED_ROOT_SNARL",
    "CHAIN", "DISTANCED_CHAIN", "MULTICOMPONENT_CHAIN",
    "CHILDREN",
    "OVERSIZED_ROOT_SNARL"
};
const vector<std::string> SnarlDistanceIndex::connectivity_t_as_string = {
    "START_START", "START_END", "START_TIP", 
    "END_START", "END_END", "END_TIP", 
    "TIP_START", "TIP_END", "TIP_TIP"
};
const vector<std::string> SnarlDistanceIndex::net_handle_record_t_string = {
    "ROOT_HANDLE", "NODE_HANDLE", "SNARL_HANDLE", 
    "CHAIN_HANDLE", "SENTINEL_HANDLE"
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//Constructor
SnarlDistanceIndex::SnarlDistanceIndex() {
    snarl_tree_records.construct(get_prefix());
    snarl_tree_records->width(64);
}
SnarlDistanceIndex::~SnarlDistanceIndex() {
}


/*Temporary distance index for constructing the index from the graph
 */
SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryDistanceIndex(){}

SnarlDistanceIndex::TemporaryDistanceIndex::~TemporaryDistanceIndex(){}

string SnarlDistanceIndex::TemporaryDistanceIndex::structure_start_end_as_string(temp_record_ref_t index) const {
    try {
        if (index.first == TEMP_NODE) {
            const TemporaryNodeRecord& temp_node_record = get_node(index);
            assert(index.second == temp_node_record.node_id);
            return "node " + std::to_string(temp_node_record.node_id);
        } else if (index.first == TEMP_SNARL) {
            const TemporarySnarlRecord& temp_snarl_record = get_snarl(index);
            return "snarl " + std::to_string(temp_snarl_record.start_node_id) 
                    + (temp_snarl_record.start_node_rev ? " rev" : " fd") 
                    + " -> " + std::to_string(temp_snarl_record.end_node_id) 
                    + (temp_snarl_record.end_node_rev ? " rev" : " fd");
        } else if (index.first == TEMP_CHAIN) {
            const TemporaryChainRecord& temp_chain_record = get_chain(index);
            return "chain " + std::to_string(temp_chain_record.start_node_id) 
                    + (temp_chain_record.start_node_rev ? " rev" : " fd") 
                    + " -> "  + std::to_string(temp_chain_record.end_node_id) 
                    + (temp_chain_record.end_node_rev ? " rev" : " fd");
        } else if (index.first == TEMP_ROOT) {
            return (string) "root";
        } else {
            return (string)"???" + std::to_string(index.first) + "???";
        }
    } catch (std::out_of_range& e) {
        throw std::out_of_range("Unable to look up (" + std::to_string(index.first) + ", " + std::to_string(index.second) + ") in temporary distance index due to out of range error: " + e.what());
    }
}
//The max record length of this chain
size_t SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryChainRecord::get_max_record_length(bool include_distances) const {
    if (is_trivial) {
        return NODE_RECORD_SIZE;
    } else { 
        //Count how many trivial snarls are in this chain
        size_t trivial_snarl_count = 0;
        size_t last_node_count = 0; // How many nodes have we seen in a row?
        size_t nontrivial_snarl_count = 0;
        size_t total_node_count = 0;
        for (const temp_record_ref_t& child : children) {
            if (child.first == TEMP_NODE) {
                if (total_node_count==0 || child != children.front()) {
                    //If this is the last node in the chain, don't do anything
                    total_node_count++;
                    if (last_node_count == 0 || last_node_count == MAX_TRIVIAL_SNARL_NODE_COUNT) {
                        //New trivial snarl
                        trivial_snarl_count ++;
                        last_node_count = 1;
                    } else {
                        last_node_count ++;
                    }
                }
            } else {
                //IF this is a snarl
                last_node_count = 0;
                nontrivial_snarl_count++;
            }
        }
        //The size of the chain record + the size of all the trivial snarls in the chain
        // + the size of the nodes in the trivial snarls + the sizes of all
        // the snarls in the chain
        return CHAIN_RECORD_SIZE + ((include_distances ? DISTANCED_TRIVIAL_SNARL_RECORD_SIZE : DISTANCELESS_TRIVIAL_SNARL_RECORD_SIZE) *trivial_snarl_count) 
                                 + (total_node_count * (include_distances ? 2 : 1)) 
                                 + ((trivial_snarl_count + nontrivial_snarl_count) * 2) - 1;
    }
}
//The max record length of the root
size_t SnarlDistanceIndex::TemporaryDistanceIndex::get_max_record_length() const {
    return ROOT_RECORD_SIZE + root_structure_count + (max_node_id-min_node_id+1) * 2 + max_index_size;
}

//The max record length of this snarl
size_t SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord::get_max_record_length() const {
    if (is_trivial) {
        return 0;
    } else if (is_simple) {
        return SimpleSnarlRecord::record_size(node_count, include_distances); 
    } else {
         if (parent.first == TEMP_ROOT) {
             if (!hub_labels.empty()) {
                 return SnarlRecord::record_size(encode_root_snarl(true, true), node_count, hub_labels.size());
             }
             //TODO: Why is node_count being added?
             return SnarlRecord::record_size(encode_root_snarl(include_distances), node_count, 0) + node_count;
         } else if (!(hub_labels.empty())) {
            return SnarlRecord::record_size(encode_nonroot_nonsimple_snarl(true, is_regular, true), node_count, hub_labels.size());
         } else {
            //TODO: Why is node_count being added?
            return SnarlRecord::record_size(encode_nonroot_nonsimple_snarl(include_distances, is_regular, false), node_count, 0) + node_count;
         }
    }
}

void SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord::allocate_distances() {
    size_t side_count = node_count * 2;
    distances.assign(side_count * (side_count + 1) / 2, std::numeric_limits<size_t>::max());
}

size_t SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord::distance_offset(size_t rank1, bool right_side1,
        size_t rank2, bool right_side2) const {
    // Same side numbering as the final index: internal children of a non-root snarl start at rank 2
    size_t first_rank = is_root_snarl ? 0 : 2;
    size_t side1 = (rank1 - first_rank) * 2 + right_side1;
    size_t side2 = (rank2 - first_rank) * 2 + right_side2;
    if (side1 > side2) {
        std::swap(side1, side2);
    }
    // Row side1 of the upper triangle starts after rows 0..side1-1, which hold side_count, side_count-1, ... entries
    size_t side_count = node_count * 2;
    return side1 * (2 * side_count - side1 + 1) / 2 + (side2 - side1);
}

void SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord::for_each_distance(
        const std::function<void(size_t, bool, size_t, bool, size_t)>& iteratee) const {
    size_t first_rank = is_root_snarl ? 0 : 2;
    size_t side_count = node_count * 2;
    size_t offset = 0;
    for (size_t side1 = 0; side1 < side_count && offset < distances.size(); side1++) {
        for (size_t side2 = side1; side2 < side_count; side2++, offset++) {
            if (distances[offset] != std::numeric_limits<size_t>::max()) {
                iteratee(side1 / 2 + first_rank, side1 % 2, side2 / 2 + first_rank, side2 % 2, distances[offset]);
            }
        }
    }
}





///////////////////////////////////////////////////////////////////////////////////////////////////
//Implement the SnarlDecomposition's functions for moving around the snarl tree
//


net_handle_t SnarlDistanceIndex::get_root() const {
    // The root is the first thing in the index, the traversal is tip to tip
    return get_net_handle_from_values(0, TIP_TIP, ROOT_HANDLE);
}

bool SnarlDistanceIndex::is_root(const net_handle_t& net) const {
    return get_handle_type(net) == ROOT_HANDLE;
}

bool SnarlDistanceIndex::is_root_snarl(const net_handle_t& net) const {

    return get_handle_type(net) == ROOT_HANDLE && get_record_offset(net) != 0;
}


bool SnarlDistanceIndex::is_snarl(const net_handle_t& net) const {
    return get_handle_type(net) == SNARL_HANDLE;
}

bool SnarlDistanceIndex::is_oversized_snarl(const net_handle_t& net) const {
    return is_oversized_snarl(SnarlTreeRecord(net, &snarl_tree_records).get_record_type());
}

bool SnarlDistanceIndex::is_dag(const net_handle_t& snarl) const {
    record_t record_type = SnarlTreeRecord(snarl, &snarl_tree_records).get_record_type();
    if ( is_nontrivial_snarl(record_type) && !has_distances(record_type) ) {
        //If this is a snarl but didn't store distances
        cerr << "warning: checking if a snarl is a dag in an index without distances. Returning true" << endl;
        return true; 
    } else if (is_nonsimple_snarl(record_type)) {
        //If this is any kind of non-simple snarl
        //(We already ruled out not having distances)

        if (!is_root_snarl(record_type)) {
            //If there were boundary nodes, check for loops on the bounds

            //The bounds of the snarl facing in
            net_handle_t snarl_start = get_bound(snarl, false, true);
            net_handle_t snarl_end = get_bound(snarl, true, true);

            //If there are loops on the bounds, then this is not a dag
            if (distance_in_parent(snarl, snarl_start, snarl_start) != std::numeric_limits<size_t>::max() ||
                distance_in_parent(snarl, snarl_end, snarl_end) != std::numeric_limits<size_t>::max()){
                return false;
            }
        }

        //Otherwise, go through each child of the snarl
        //This returns false if it stopped early - if the iterator returned false - if it had non-dag edges
        return for_each_child(snarl, [&] (const net_handle_t& child) {
            //Check if there are any self-loops on the child
            if (distance_in_parent(snarl, child, child) != std::numeric_limits<size_t>::max() ||
                distance_in_parent(snarl, child, flip(child)) != std::numeric_limits<size_t>::max() ||
                distance_in_parent(snarl, flip(child), child) != std::numeric_limits<size_t>::max()){

                //If there are self-loops
                return false;
            } else {
                return true;
            }
        });
    } else {
        //If this is a simple snarl or not a snarl, then return false
        if (!is_snarl(snarl)) {
            cerr << "Warning: checking if a non-snarl is a dag" << endl;
        }
        return true;
    }
}

bool SnarlDistanceIndex::is_simple_snarl(const net_handle_t& net) const {
    return get_handle_type(net) == SNARL_HANDLE && get_node_record_offset(net) == 1;
}

bool SnarlDistanceIndex::is_regular_snarl(const net_handle_t& net) const {
    record_t record_type = SnarlTreeRecord(net, &snarl_tree_records).get_record_type();
    return is_regular_snarl(record_type);
}

bool SnarlDistanceIndex::is_chain(const net_handle_t& net) const {
    return get_handle_type(net) == CHAIN_HANDLE;
}
bool SnarlDistanceIndex::is_multicomponent_chain(const net_handle_t& net) const {
    return get_handle_type(net) == CHAIN_HANDLE 
        && SnarlTreeRecord(net, &snarl_tree_records).get_record_type() == MULTICOMPONENT_CHAIN;
}
bool SnarlDistanceIndex::is_looping_chain(const net_handle_t& net) const {
    if (!is_chain(net) || is_trivial_chain(net)) {
        return false;
    }
    ChainRecord chain_record(net, &snarl_tree_records);
    return chain_record.get_start_id() == chain_record.get_end_id();
}
bool SnarlDistanceIndex::is_ordered_in_chain(const net_handle_t& child1, const net_handle_t& child2) const {
    size_t rank1 = get_record_offset(child1) + get_node_record_offset(child1);
    size_t rank2 = get_record_offset(child2) + get_node_record_offset(child2);

    return  rank1 < rank2;
}


bool SnarlDistanceIndex::is_trivial_chain(const net_handle_t& net) const {
    bool handle_is_chain =get_handle_type(net) == CHAIN_HANDLE; 
    bool record_is_node = SnarlTreeRecord(net, &snarl_tree_records).get_record_handle_type() == NODE_HANDLE;
    bool record_is_simple_snarl = is_simple_snarl(SnarlTreeRecord(net, &snarl_tree_records).get_record_type());
    bool handle_has_node_offset = get_node_record_offset(net) >= 2;
    
    return handle_is_chain && (record_is_node
                || (record_is_simple_snarl && handle_has_node_offset));
}

bool SnarlDistanceIndex::is_node(const net_handle_t& net) const {
    return get_handle_type(net) == NODE_HANDLE;
}
bool SnarlDistanceIndex::is_sentinel(const net_handle_t& net) const {
    return get_handle_type(net) == SENTINEL_HANDLE;

}

net_handle_t SnarlDistanceIndex::get_net(const handle_t& handle, const handlegraph::HandleGraph* graph) const{
    net_handle_t net = get_node_net_handle(graph->get_id(handle));
    return graph->get_is_reverse(handle) ? flip(net) : net;

}
handle_t SnarlDistanceIndex::get_handle(const net_handle_t& net, const handlegraph::HandleGraph* graph) const{
    if (get_handle_type(net) == SENTINEL_HANDLE) {
        // We don't know if this is a trivial or nontrivial snarl, so we need
        // to access it with the base class.
        // TODO: Make a base class for any kind of snarl.
        SnarlTreeRecord snarl_record(net, &snarl_tree_records);
        if (starts_at(net) == START) {
            return graph->get_handle(snarl_record.get_start_id(), 
                       ends_at(net) == START ? !snarl_record.get_start_orientation()   //Going out
                                             : snarl_record.get_start_orientation());  //Going in
        } else {
            assert (starts_at(net) == END);
            return graph->get_handle(snarl_record.get_end_id(), 
                       ends_at(net) == END ? snarl_record.get_end_orientation()   //Going out
                                           : !snarl_record.get_end_orientation());  //Going in
        }

    } else if (get_handle_type(net) == NODE_HANDLE || get_handle_type(net) == CHAIN_HANDLE ) {
        return graph->get_handle(node_id(net), ends_at(net) == END ? false : true);
    } else {
        throw runtime_error("error: trying to get a handle from a snarl, chain, or root");
    }
}

net_handle_t SnarlDistanceIndex::get_parent(const net_handle_t& child) const {

    connectivity_t child_connectivity = get_connectivity(child);

    //If the child is the sentinel of a snarl, just return the snarl
    if (get_handle_type(child) == SENTINEL_HANDLE) {
        connectivity_t parent_connectivity = START_END;
        if (child_connectivity == START_END || child_connectivity == END_START) {
            //If the child sentinel is pointing through the snarl, then the snarl should point
            //in the same direction
            parent_connectivity = child_connectivity;
        } else if (child_connectivity == START_START) {
            //If the sentinel is the start node pointing out, then the parent should point out the start
            parent_connectivity = END_START;
        } else if (child_connectivity == END_END) {
            //If the sentinel is the end node pointing out, then the parent should point out the end
            parent_connectivity = START_END;
        }
        return get_net_handle_from_values(get_record_offset(child), parent_connectivity, 
                                          SNARL_HANDLE, get_node_record_offset(child)); 
    } else if (get_handle_type(child) == ROOT_HANDLE) {
        throw runtime_error("error: trying to find the parent of the root");
    } else if (is_simple_snarl(get_record_type(snarl_tree_records->at(get_record_offset(child))))) {

        //If this is a simple snarl and a node or chain, then the parent offset doesn't change
        if (get_handle_type(child) == NODE_HANDLE) {
            // If this is a node, then return it as a chain
            return get_net_handle_from_values(get_record_offset(child), child_connectivity, 
                                              CHAIN_HANDLE, get_node_record_offset(child));
        } else if (get_handle_type(child) == CHAIN_HANDLE) {
            //If this is a chain, then return the same thing as a snarl
            return get_net_handle_from_values(get_record_offset(child), START_END, SNARL_HANDLE, 1);
        }
    }

    //Otherwise, we need to move up one level in the snarl tree

    SnarlTreeRecord child_record(child, &snarl_tree_records);
    //Get the pointer to the parent to find its type
    size_t parent_pointer = child_record.get_parent_record_offset();
    SnarlTreeRecord parent_record(parent_pointer, &snarl_tree_records);
    net_handle_record_t parent_type = parent_record.get_record_handle_type();

    
    //The connectivity of the parent defaults to start-end
    connectivity_t parent_connectivity = START_END;
    //If the parent is going to be a chain, then it has the connectivity of the child, flipped if the child is 
    //reversed in the parent
    if ((get_handle_type(child) == NODE_HANDLE || get_handle_type(child) == SNARL_HANDLE) 
            && (child_connectivity == START_END || child_connectivity == END_START)) {
        if (is_reversed_in_parent(child)) {
            parent_connectivity = child_connectivity == START_END ? END_START : START_END;

        } else {
            parent_connectivity = child_connectivity;
        }
    } else if (parent_type == ROOT_HANDLE) {
        parent_connectivity = TIP_TIP;
    }
    if (get_handle_type(child) == NODE_HANDLE && parent_type != CHAIN_HANDLE) {
        //If this is a node and it's parent is not a chain, we want to pretend that its 
        //parent is a chain version of the child
        net_handle_t projected = get_net_handle_from_values(get_record_offset(child), child_connectivity, CHAIN_HANDLE, get_node_record_offset(child));

        return projected;
    } 

    return get_net_handle(parent_pointer, parent_connectivity);
}

net_handle_t SnarlDistanceIndex::get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const {
    if (get_handle_type(snarl) == CHAIN_HANDLE) {
        // This could be a real chain, a node looking like a chain, or a simple
        // snarl record looking like a chain (maybe because the node it was
        // looking like needs to look like a chain now). ChainRecord promises
        // to know how to interpret all of them.
        // TODO: the concepts involved in things looking like other things
        // should be documented somewhere.
        ChainRecord chain_record(snarl, &snarl_tree_records);
        size_t offset;
        size_t node_offset;
        bool is_looping_chain=false;
        bool rev_in_parent;
        if (is_trivial_chain(snarl)) {
            offset = get_record_offset(snarl);
            rev_in_parent = false;
            node_offset = get_node_record_offset(snarl);
        } else {
            if (get_end) {
                tuple<size_t, bool, size_t> last_child = chain_record.get_last_child_offset();
                is_looping_chain = chain_record.get_start_id() == chain_record.get_end_id();
                offset = is_looping_chain ? chain_record.get_first_node_offset() : std::get<0>(last_child);
                node_offset = is_looping_chain ? 0 : std::get<2>(last_child); 

            } else {
                offset = chain_record.get_first_node_offset(); 
                node_offset = 0;
                
            }
            rev_in_parent = TrivialSnarlRecord(offset, &snarl_tree_records).get_is_reversed_in_parent(node_offset);
        }
        if (get_end) {
            rev_in_parent = !rev_in_parent;
        }
        if (!face_in){
            rev_in_parent = !rev_in_parent;
        }
        connectivity_t connectivity = rev_in_parent ? END_START : START_END;

        if (is_looping_chain && get_end) {
            //If this is a looping chain and we're getting the end, then the traversal of the node will be the 
            //end endpoint repeated, instead of the actual traversal
            connectivity = endpoints_to_connectivity(get_end_endpoint(connectivity), get_end_endpoint(connectivity));
        }
        return get_net_handle_from_values(offset, connectivity,  NODE_HANDLE, node_offset);

    } else if (is_root(snarl)) {
        throw runtime_error("Trying to get the bounds of a root");
    } else {
        assert(get_handle_type(snarl) == SNARL_HANDLE);
        endpoint_t start = get_end ? END : START;
        endpoint_t end = face_in ? (start == END ? START : END) : start;
        return get_net_handle_from_values(get_record_offset(snarl), endpoints_to_connectivity(start, end), SENTINEL_HANDLE, get_node_record_offset(snarl));
    }
}

net_handle_t SnarlDistanceIndex::get_node_from_sentinel(const net_handle_t& sentinel) const {
    assert(is_sentinel(sentinel));
    if (is_root(get_parent(get_parent(sentinel)))) {
        throw runtime_error("Trying to get the bounds of a root");
    }
    net_handle_t snarl = get_parent(sentinel);
    net_handle_t parent_chain = get_parent(snarl);
    assert(is_chain(parent_chain));
    //The snarl must be in a chain, so the boundary nodes will be the next things on the chain 
    ChainRecord chain_record(parent_chain, &snarl_tree_records);

    net_handle_t next_handle = chain_record.get_next_child(snarl, starts_at(sentinel) == START);

    //If the sentinel is facing into the snarl, then flip the node after getting it from the chain
    return ends_at(sentinel) != starts_at(sentinel) ? flip(next_handle) : next_handle ;


}

net_handle_t SnarlDistanceIndex::flip(const net_handle_t& net) const {
    connectivity_t old_connectivity = get_connectivity(net);
    connectivity_t new_connectivity;
    if (is_sentinel(net)) {
        //If this is the sentinel of a snarl, then the first endpoint indicates the node
        // and the second indicates the orientation, so only flip the second
        new_connectivity = endpoints_to_connectivity(get_start_endpoint(old_connectivity),
                get_end_endpoint(old_connectivity) == START ? END : START);
    } else {
        new_connectivity = endpoints_to_connectivity(get_end_endpoint(old_connectivity), 
                                                                get_start_endpoint(old_connectivity));
    }
    return get_net_handle_from_values(get_record_offset(net), new_connectivity, get_handle_type(net), get_node_record_offset(net));
}

net_handle_t SnarlDistanceIndex::canonical(const net_handle_t& net) const {
    SnarlTreeRecord record(net, &snarl_tree_records);
    record_t type = record.get_record_type();
    if (is_any_root(type)) {
        return get_root();
    }


    connectivity_t connectivity;
    if (record.is_start_end_connected()) {
        connectivity = START_END;
    } else if (record.is_start_tip_connected()) {
        connectivity = START_TIP;
    } else if (record.is_end_tip_connected()) {
        connectivity = END_TIP;
    } else if (record.is_start_start_connected()) {
        connectivity = START_START;
    } else if (record.is_end_end_connected()) {
        connectivity = END_END;
    } else if (record.is_tip_tip_connected()) {
        connectivity = TIP_TIP;
    } else {
        connectivity = START_END; 
    }
    return get_net_handle_from_values(get_record_offset(net), connectivity, get_handle_type(net), get_node_record_offset(net));
}

net_handle_t SnarlDistanceIndex::start_end_traversal_of(const net_handle_t& net) const {
    return get_net_handle_from_values(get_record_offset(net), START_END, get_handle_type(net), get_node_record_offset(net));
}

SnarlDecomposition::endpoint_t SnarlDistanceIndex::starts_at(const net_handle_t& traversal) const {
    return get_start_endpoint(get_connectivity(traversal));

}
SnarlDecomposition::endpoint_t SnarlDistanceIndex::ends_at(const net_handle_t& traversal) const {
    return get_end_endpoint( get_connectivity(traversal));
}

size_t SnarlDistanceIndex::get_rank_in_parent(const net_handle_t& net) const {
    size_t tag = snarl_tree_records->at(get_record_offset(net));
    if (is_trivial_snarl(get_record_type(tag))) {
        return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_rank_in_parent(get_node_record_offset(net));
    } else if (is_simple_snarl(get_record_type(tag))) {
        if (is_snarl(net)) {
            return get_record_offset(net);
        } else {
            return get_node_record_offset(net);
        }
    } else {
        return SnarlTreeRecord(net, &snarl_tree_records).get_rank_in_parent();
    }
}

size_t SnarlDistanceIndex::connected_component_count() const {
    return RootRecord (get_root(), &snarl_tree_records).get_connected_component_count();
}

net_handle_t SnarlDistanceIndex::get_snarl_child_from_rank(const net_handle_t& snarl, const size_t& rank) const {
    if (rank == 0) {
        return get_bound(snarl, false, true);
    } else if (rank == 1) {
        return get_bound(snarl, true, true);
    } else if (is_simple_snarl(snarl) ){
        return SimpleSnarlRecord(snarl, &snarl_tree_records).get_child_from_rank(rank);
    } else {
        //Ranks for children of snarls start from 2 since 0 and 1 are reserved for the bounds
        if (rank == 0) {
            //The boundary node. This technically shouldn't be here because boundary nodes aren't children
            //of the snarl but I'll leave it anyway because I will probably use it
            //Returns the bound facing in
            return get_bound(snarl, false, true);
        } else if ( rank == 1) {
            //The end bound
            return get_bound(snarl, true, true);
        } else {
            //TODO: Since I changed the ranks of snarls to be in a topological order, this is now slow
            net_handle_t child_handle;
            for_each_child(snarl, [&](const net_handle_t& child) {
                if (get_rank_in_parent(child) == rank) {
                    child_handle = child;
                    return false;
                }
                return true;
            });
            return child_handle;
        }
    }
}

bool SnarlDistanceIndex::has_distances(const net_handle_t& net) const {
    return has_distances(SnarlTreeRecord(net, &snarl_tree_records).get_record_type()); 
}

bool SnarlDistanceIndex::has_distances() const {
    RootRecord root_record (get_root(), &snarl_tree_records);
    return has_distances(get_node_net_handle(root_record.get_min_node_id())); 
}

size_t SnarlDistanceIndex::get_snarl_child_count(const net_handle_t& net) const {
    record_t specific_type = SnarlTreeRecord(net, &snarl_tree_records).get_record_type();
    if (is_simple_snarl(specific_type)) {
        return SimpleSnarlRecord(net, &snarl_tree_records).get_node_count();
    } else {
        return SnarlRecord(net, &snarl_tree_records).get_node_count();
    }
}

bool SnarlDistanceIndex::for_each_child_impl(const net_handle_t& traversal, const std::function<bool(const net_handle_t&)>& iteratee) const {
    //What is this according to the snarl tree
    net_handle_record_t record_type = SnarlTreeRecord(traversal, &snarl_tree_records).get_record_handle_type();
    //What is this according to the handle 
    //(could be a trivial chain but actually a node according to the snarl tree)
    net_handle_record_t handle_type = get_handle_type(traversal);
    if (record_type == ROOT_HANDLE) {
        RootRecord root_record(get_root(), &snarl_tree_records);
        return root_record.for_each_child(iteratee);
    } else if (is_simple_snarl(SnarlTreeRecord(traversal, &snarl_tree_records).get_record_type())) {
        //If this is a simple snarl then it is a bit different
        if (handle_type == CHAIN_HANDLE) {
            //If the handle thinks it's a chain, then it is a trivial chain in the snarl so we do
            //iteratee for the node
            return iteratee(get_net_handle_from_values(get_record_offset(traversal), get_connectivity(traversal), 
                                           NODE_HANDLE, get_node_record_offset(traversal)));
        } else if (handle_type == SNARL_HANDLE) {
            return SimpleSnarlRecord(traversal, &snarl_tree_records).for_each_child(iteratee);
        } else { 
            throw runtime_error("error: Looking for children of a node or sentinel in a simple snarl");
        }
    } else if (record_type == SNARL_HANDLE) {
        // This could be a simple or non-simple snarl
        record_t specific_type = SnarlTreeRecord(traversal, &snarl_tree_records).get_record_type();
        if (is_simple_snarl(specific_type)) {
            SimpleSnarlRecord snarl_record(traversal, &snarl_tree_records);
            return snarl_record.for_each_child(iteratee);
        } else {
            SnarlRecord snarl_record(traversal, &snarl_tree_records);
            return snarl_record.for_each_child(iteratee);
        }
    } else if (record_type == CHAIN_HANDLE) {
        ChainRecord chain_record(traversal, &snarl_tree_records);
        return chain_record.for_each_child(iteratee);
    } else  if (record_type == NODE_HANDLE && handle_type == CHAIN_HANDLE) {
        //This is actually a node but we're pretending it's a chain
        return iteratee(get_net_handle_from_values(get_record_offset(traversal), get_connectivity(traversal), NODE_HANDLE));
    } else {
        throw runtime_error("error: Looking for children of a node or sentinel");
    }
   
}

bool SnarlDistanceIndex::for_each_traversal_impl(const net_handle_t& item, const std::function<bool(const net_handle_t&)>& iteratee) const {
    if (get_handle_type(item) == SENTINEL_HANDLE) {
        if (!iteratee(get_net_handle_from_values(get_record_offset(item), START_END, get_handle_type(item), get_node_record_offset(item)))) {
            return false;
        }
        if (!iteratee(get_net_handle_from_values(get_record_offset(item), END_START, get_handle_type(item), get_node_record_offset(item)))) {
            return false;
        }
        return true;
    }
    SnarlTreeRecord record(item, &snarl_tree_records);
    for ( size_t type = 1 ; type <= 9 ; type ++ ){
        connectivity_t connectivity = static_cast<connectivity_t>(type);
        if (record.has_connectivity(connectivity)) {
            if (!iteratee(get_net_handle_from_values(get_record_offset(item), connectivity, get_handle_type(item), get_node_record_offset(item)))) {
                return false;
            }
        }
    }
    return true;
}

bool SnarlDistanceIndex::follow_net_edges_impl(const net_handle_t& here, const handlegraph::HandleGraph* graph, bool go_left, const std::function<bool(const net_handle_t&)>& iteratee) const {

    SnarlTreeRecord this_record(here, &snarl_tree_records);
    SnarlTreeRecord parent_record (get_parent(here), &snarl_tree_records);

    if (parent_record.get_record_handle_type() == ROOT_HANDLE &&
        !is_root_snarl(parent_record.get_record_type())) {
        // TODO: should we check for ROOT record type here?
        //If this is a root-level chain or node
        if ((go_left && starts_at(here) == END) || (!go_left && ends_at(here) == END)) {
            //Follow edges leaving the root structure at the end
            if (this_record.is_externally_start_end_connected()) {
                //Follow edge from end to start
                if (!iteratee(get_net_handle_from_values(get_record_offset(here), START_END, get_handle_type(here), get_node_record_offset(here)))) {
                    return false;
                }
            }
            if (this_record.is_externally_end_end_connected()) {
                //Follow edge from end back to the end
                if (!iteratee(get_net_handle_from_values(get_record_offset(here), END_START, get_handle_type(here), get_node_record_offset(here)))) {
                    return false;
                }
            }
        } else {
            //Follow edges leaving the root structure at the end
            if (this_record.is_externally_start_end_connected()) {
                //Follow edge from start to end
                if (!iteratee(get_net_handle_from_values(get_record_offset(here), END_START, get_handle_type(here), get_node_record_offset(here)))) {
                    return false;
                }
            }
            if (this_record.is_externally_end_end_connected()) {
                //Follow edge from the start back to the start
                if (!iteratee(get_net_handle_from_values(get_record_offset(here), START_END, get_handle_type(here), get_node_record_offset(here)))) {
                    return false;
                }
            }

        }
        return true;

    } else if (get_handle_type(here) == CHAIN_HANDLE || get_handle_type(here) == SENTINEL_HANDLE) {
        //If this is a chain (or a node pretending to be a chain) and it is the child of a snarl
        //Or if it is the sentinel of a snarl, then we walk through edges in the snarl
        //It can either run into another chain (or node) or the boundary node
        bool is_root_snarl = SnarlDistanceIndex::is_root_snarl(parent_record.get_record_type());


        //Get the graph handle for the end node of whatever this is, pointing in the right direction
        handle_t graph_handle;
        if (get_handle_type(here) == SENTINEL_HANDLE) {
            if ((get_connectivity(here) == START_END && !go_left) ||
                (get_connectivity(here) == START_START && go_left)) {
                //If it's going into the snarl from the start node
                graph_handle = graph->get_handle(parent_record.get_start_id(), parent_record.get_start_orientation());
            } else if ((get_connectivity(here) == END_START && !go_left) ||
                       (get_connectivity(here) == END_END && go_left)) {
                //If it's going into the snarl from the end node
                graph_handle = graph->get_handle(parent_record.get_end_id(), !parent_record.get_end_orientation());
            } else {
                //This is facing out, so don't do anything 
                return true;
            }
        } else if (get_handle_type(here) == NODE_HANDLE ||is_trivial_chain(here)) {
            nid_t id = node_id(here);
            graph_handle = graph->get_handle(id, go_left ? starts_at(here) == START : ends_at(here) == START);
        } else {
            //Get the bound of the snarl or chain. If the traversal ends (or, if go_left is true, starts) at the
            //start, then go left from here, otherwise go right (forwards)
            graph_handle = get_handle(get_bound(here, go_left ? starts_at(here) != START : ends_at(here) != START, false), graph);
        }

        // Follow edges in the graph. Return the return value of graph follow_edges
        return graph->follow_edges(graph_handle, false, [&](const handle_t& h) {

            if (!is_root_snarl && graph->get_id(h) == parent_record.get_start_id()) {
                //If this is the start boundary node of the parent snarl, then do this on the sentinel
                return iteratee(get_bound(get_parent(here), false, false));
            } else if (!is_root_snarl && graph->get_id(h) == parent_record.get_end_id()) {
                return iteratee(get_bound(get_parent(here), true, false));
            } else {
                //It is either another chain or a node, but the node needs to pretend to be a chain
                net_handle_t node_net_handle = get_node_net_handle(graph->get_id(h));
                if (graph->get_is_reverse(h)) {
                    node_net_handle = flip(node_net_handle);
                }

                if (is_node(get_record_type(snarl_tree_records->at(get_record_offset(node_net_handle))))) {
                    //If this is a node make a net_handle_t of a node pretending to be a chain
                    net_handle_t next_net = get_net_handle_from_values(get_record_offset(node_net_handle), 
                                                           graph->get_is_reverse(h) ? END_START : START_END, 
                                                           CHAIN_HANDLE);
                   return iteratee(next_net);
                } else if (is_simple_snarl(get_record_type(snarl_tree_records->at(get_record_offset(node_net_handle))))) {
                    //If the node is a node in a simple snarl
                    net_handle_t next_net = get_net_handle_from_values(get_record_offset(node_net_handle), 
                                                           graph->get_is_reverse(h) ? END_START : START_END, 
                                                           CHAIN_HANDLE,
                                                           get_node_record_offset(node_net_handle));
                   return iteratee(next_net);
                } else {
                    //next_node_record is also the start of a chain
                    bool rev = node_net_handle != get_bound(get_parent(node_net_handle), false, true);
                    net_handle_t next_net = get_net_handle_from_values(get_record_offset(get_parent(node_net_handle)), 
                                                           rev ? END_START : START_END, 
                                                           CHAIN_HANDLE);
                   return iteratee(next_net);
                }
            }
            return false;
        });

        
    } else if (get_handle_type(here) == SNARL_HANDLE || get_handle_type(here) == NODE_HANDLE) {
        if (get_record_offset(here) == parent_record.record_offset) {
            //If this is a trivial chain
            return true;
        }
        //If this is a snarl or node, then it is the component of a (possibly pretend) chain
        ChainRecord parent_chain(this_record.get_parent_record_offset(), &snarl_tree_records);
        bool go_left_in_chain = go_left ? starts_at(here) == START : ends_at(here) == START;
        bool is_rev = is_node(here) ? TrivialSnarlRecord(get_record_offset(here), &snarl_tree_records).get_is_reversed_in_parent(get_node_record_offset(here))
                                    : false;
        if (is_rev) {
            go_left_in_chain = !go_left_in_chain;
        }
        net_handle_t next_net = parent_chain.get_next_child(here, go_left_in_chain);
        if (next_net == here ) {
            //If this is the end of the chain
            return true;
        }
        return iteratee(next_net);
        
    }
    return true;
}


net_handle_t SnarlDistanceIndex::get_parent_traversal(const net_handle_t& traversal_start, const net_handle_t& traversal_end) const {
    
    net_handle_record_t start_handle_type = get_handle_type(traversal_start);
    net_handle_record_t end_handle_type = get_handle_type(traversal_end);

    if (start_handle_type == SENTINEL_HANDLE) {
        //these are the sentinels of a snarl
        assert(end_handle_type == SENTINEL_HANDLE);
        endpoint_t start_endpoint = get_start_endpoint(get_connectivity(traversal_start));
        endpoint_t end_endpoint = get_start_endpoint(get_connectivity(traversal_end));
        return get_net_handle_from_values(get_record_offset(get_parent(traversal_start)), 
                              endpoints_to_connectivity(start_endpoint, end_endpoint),
                              SNARL_HANDLE, get_node_record_offset(traversal_start));
    } else {
        //These are the endpoints or tips in a chain
        SnarlTreeRecord start_record = get_snarl_tree_record(traversal_start);
        SnarlTreeRecord end_record = get_snarl_tree_record(traversal_end);
        SnarlTreeRecord parent_record (start_record.get_parent_record_offset(), &snarl_tree_records);

        //Figure out what the start and end of the traversal are
        endpoint_t start_endpoint;
        if (start_handle_type == NODE_HANDLE && 
            node_id(traversal_start) == parent_record.get_start_id() &&
            (get_start_endpoint(traversal_start) == START && !parent_record.get_start_orientation() ||
             get_start_endpoint(traversal_start) == END && parent_record.get_start_orientation()) ){
            //If traversal_start is a node and is also the start node oriented into the parent
            start_endpoint = START;
    
        } else if (start_handle_type == NODE_HANDLE && 
            node_id(traversal_start) == parent_record.get_end_id() &&
            (get_start_endpoint(traversal_start) == START && parent_record.get_end_orientation() ||
             get_start_endpoint(traversal_start) == END && !parent_record.get_end_orientation()) ){
            //If traversal_start is a node and also the end node and oriented going into the parent
            start_endpoint = END;
    
        } else if (start_handle_type == SNARL_HANDLE) {
            start_endpoint = TIP;
        } else {
            throw runtime_error("error: trying to get an invalid traversal of a chain");
        }
    
        endpoint_t end_endpoint;
        if (end_handle_type == NODE_HANDLE && 
            node_id(traversal_end) == parent_record.get_start_id() &&
            (get_start_endpoint(traversal_end) == START && parent_record.get_start_orientation() ||
             get_start_endpoint(traversal_end) == END && !parent_record.get_start_orientation())){
            //If traversal_end is a node and also the start node oriented out of the parent
            end_endpoint = START;
        } else if (end_handle_type == NODE_HANDLE && 
            node_id(traversal_end) == parent_record.get_end_id() &&
            (get_start_endpoint(traversal_end) == START && !parent_record.get_end_orientation() ||
             get_start_endpoint(traversal_end) == END && parent_record.get_end_orientation()) ){
            //If traversal_end is a node and also the end node oriented out of the parent
            end_endpoint = END;
        } else if (end_handle_type == SNARL_HANDLE) {
            end_endpoint = TIP;
        } else {
            throw runtime_error("error: trying to get an invalid traversal of a chain");
        }

        return get_net_handle_from_values(parent_record.record_offset, 
                              endpoints_to_connectivity(start_endpoint, end_endpoint), 
                              CHAIN_HANDLE);
    }


}


//Serialize and deserialize from TriviallySerializable will just call save and load

void SnarlDistanceIndex::dissociate() {
    snarl_tree_records.dissociate();
    snarl_tree_records.reset();
}
void SnarlDistanceIndex::serialize(const std::function<void(const void*, size_t)>& iteratee)  const {
    snarl_tree_records.save(iteratee);
}
void SnarlDistanceIndex::serialize(int fd) {
    //This gets called by TriviallySerializable::serialize(filename), which doesn't
    //write the prefix so it must be written here
    snarl_tree_records.save(fd);
}
void SnarlDistanceIndex::deserialize(int fd) {
    //This gets called by TriviallySerializable::deserialize(filename), which
    //doesn't check for the prefix, so this should expect it
    snarl_tree_records.load(fd, get_prefix());
    check_version_on_load();
}

void SnarlDistanceIndex::serialize_members(std::ostream& out) const {
    //This gets called by Serializable::serialize(ostream), which writes the prefix
    //so don't write the prefix here
    snarl_tree_records.save_after_prefix(out, get_prefix());
}
void SnarlDistanceIndex::deserialize_members(std::istream& in){
    //This gets called by Serializable::deserialize(istream), which has already
    //read the prefix, so don't expect the prefix
    snarl_tree_records.load_after_prefix(in, get_prefix());
    check_version_on_load();
}

void SnarlDistanceIndex::check_version_on_load() const {
    RootRecord root_record (get_root(), &snarl_tree_records);
    if (root_record.get_version_number() != CURRENT_VERSION_NUMBER) {
        throw runtime_error("error: Trying to load a SnarlDistanceIndex which is v" +
                            std::to_string(root_record.get_version_number()) + " instead of the up-to-date v" +
                            std::to_string(CURRENT_VERSION_NUMBER) + "; try regenerating the index");
    }
}

uint32_t SnarlDistanceIndex::get_magic_number()const {
    //random number?
    return magic_number;
}
//Copied from MappedPackedGraph
std::string SnarlDistanceIndex::get_prefix() const {
    //return the magic number as a string

    //Reinterpret the magic number as a char array (with 4 chars) 
    uint32_t number = htonl(magic_number);
    const char* as_chars = reinterpret_cast<const char*> (&number);
    std::string as_string(as_chars, as_chars+4);
    return as_string;
}

void SnarlDistanceIndex::preload(bool blocking) const {
    snarl_tree_records.preload(blocking);
}

size_t SnarlDistanceIndex::distance_in_parent(const net_handle_t& parent, 
        const net_handle_t& child1, const net_handle_t& child2, const HandleGraph* graph, size_t distance_limit) const {


    //Get the orientation of the children. This only cares about the end endpoint, and assumes that things that end
    //in tips are start-end, meaning that they don't end at the start
    bool child_ends_at_start1 = ends_at(child1) == START;
    bool child_ends_at_start2 = ends_at(child2) == START;

    if (is_root(parent)) {
        //If the parent is the root, then the children must be in the same root snarl for them to be
        //connected
        size_t parent_record_offset1 = get_record_offset(get_parent(child1));
        size_t parent_record_offset2 = get_record_offset(get_parent(child2));

        if (parent_record_offset1 != parent_record_offset2) {
            //If the children are in different connected components
            return std::numeric_limits<size_t>::max();
        } else {
            record_t root_record_type = get_record_type(snarl_tree_records->at(parent_record_offset1));
            if (root_record_type != DISTANCED_ROOT_SNARL && root_record_type != OVERSIZED_ROOT_SNARL) {
                //If they are in the same connected component, but it is not a root snarl
                if (get_record_offset(child1) == get_record_offset(child2)) {
                    //If they are the same child of the root but not in a snarl, then check the external connectivity
                    if (child_ends_at_start1 && child_ends_at_start2) {
                        return has_external_connectivity(child1, START, START) ? 0 : std::numeric_limits<size_t>::max();
                    } else if (!child_ends_at_start1 && !child_ends_at_start2) {
                        return has_external_connectivity(child1, END, END) ? 0 : std::numeric_limits<size_t>::max();
                    } else if ((child_ends_at_start1 && !child_ends_at_start2) ||
                                (!child_ends_at_start1 && child_ends_at_start2)) {
                        if (has_external_connectivity(child1, START, END)) {
                            //If we can take an edge around the snarl
                            return 0;
                        } else if (has_external_connectivity(child1, START, START) && has_external_connectivity(child1, END, END)) {
                            //If we can take the loops on the two ends of the snarl, walk through the snarl
                            return minimum_length(child1);
                        }
                    }
                }
                return std::numeric_limits<size_t>::max();
            } else if (root_record_type == OVERSIZED_ROOT_SNARL) {
                //Oversized root snarls store hub labels in their record. Ranks in
                //root snarls are always >= 2 (no boundary), so the rank==0 special
                //case in the non-root HHL branch reduces to zero here.
                size_t rank1 = get_rank_in_parent(child1);
                size_t rank2 = get_rank_in_parent(child2);
                bool dir1 = !child_ends_at_start1;
                bool dir2 = !child_ends_at_start2;
                auto record_it      = snarl_tree_records->begin() + parent_record_offset1;
                auto length_data_it = record_it + SNARL_RECORD_SIZE;
                size_t from_port = bgid(rank1, !dir1, true);
                size_t to_port   = bgid(rank2,  dir2, false);
                return promote_distance<size_t>(hhl_query(length_data_it, from_port, to_port));
            } else {
                //They are in the same root snarl, so find the distance between them.
                // We know this isn't a simple snarl.
                SnarlRecord snarl_record(parent_record_offset1, &snarl_tree_records);

                return snarl_record.get_distance(get_rank_in_parent(child1), !child_ends_at_start1,
                                                 get_rank_in_parent(child2), !child_ends_at_start2);
            }
        }


    } else if (is_chain(parent)) {
        if (get_record_handle_type(get_record_type(snarl_tree_records->at(get_record_offset(parent)))) == NODE_HANDLE ||
            get_record_handle_type(get_record_type(snarl_tree_records->at(get_record_offset(parent)))) == SNARL_HANDLE) {
            return std::numeric_limits<size_t>::max();
        }
        ChainRecord chain_record(parent, &snarl_tree_records);

        //These are the values needed for calculating the distance in the chain
        //They are:
        //rank in chain, go_left, node_length, prefix_sum, forward_loop, reverse_loop, component
        size_t rank_in_chain1; size_t rank_in_chain2;
        bool go_left1; bool go_left2;
        size_t node_length1; size_t node_length2;
        size_t prefix_sum1; size_t prefix_sum2;
        size_t forward_loop1 ; size_t forward_loop2;
        size_t reverse_loop1 ; size_t reverse_loop2;
        size_t component1; size_t component2;
        size_t end_component1; size_t end_component2;

        //The node lengths of snarl boundaries, to be added to the total distance
        size_t node_lengths_to_add = 0;

        //Get all the relevant chain values
        for (size_t i = 0 ; i <=1 ; i ++ ) {
            const net_handle_t& child = i == 0 ? child1 : child2;
            size_t& rank_in_chain= i == 0 ? rank_in_chain1: rank_in_chain2;
            bool& go_left_in_chain = i == 0 ? go_left1  : go_left2;
            size_t& node_length  = i == 0 ? node_length1  : node_length2;
            size_t& prefix_sum   = i == 0 ? prefix_sum1   : prefix_sum2;
            size_t& forward_loop = i == 0 ? forward_loop1 : forward_loop2;
            size_t& reverse_loop = i == 0 ? reverse_loop1 : reverse_loop2;
            size_t& component    = i == 0 ? component1    : component2;
            size_t& end_component    = i == 0 ? end_component1    : end_component2;
            bool& child_ends_at_start = i == 0 ? child_ends_at_start1 : child_ends_at_start2;

            if (is_node(child)) {

                bool go_left = is_reversed_in_parent(child);
                if (child_ends_at_start){
                    go_left = !go_left;
                }
                rank_in_chain=get_rank_in_parent(child);
                go_left_in_chain = go_left;
                node_length = minimum_length(child);
                std::tie(prefix_sum, forward_loop, reverse_loop, component) = 
                    TrivialSnarlRecord(get_record_offset(child), &snarl_tree_records).get_chain_values(get_node_record_offset(child));
                end_component = TrivialSnarlRecord(get_record_offset(child), &snarl_tree_records).get_chain_component(get_node_record_offset(child), true);

                
            } else {
                //If this is a snarl, then find the appropriate boundary node and fill it in
                assert(is_snarl(child));
                if (child_ends_at_start) {
                    //If we're going backwards in the chain, then we want the start node
                    //And turn it into a node, since it will be a sentinel
                    net_handle_t start_bound = get_node_from_sentinel(get_bound(child, false, false));

                    size_t start_length = minimum_length(start_bound);

                        //set values for first
                    rank_in_chain=get_record_offset(start_bound) + get_node_record_offset(start_bound);
                    go_left_in_chain = true;
                    node_length = start_length;
                    std::tie(prefix_sum, forward_loop, reverse_loop, component) = 
                        TrivialSnarlRecord (get_record_offset(start_bound), &snarl_tree_records).get_chain_values(get_node_record_offset(start_bound));

                    node_lengths_to_add += start_length;
                } else {
                    //Do the same thing for the snarl end node if we're going forwards
                    net_handle_t end_bound = get_node_from_sentinel(get_bound(child, true, false));
                    size_t end_length = minimum_length(end_bound);
                    size_t end_rank = get_record_offset(end_bound) + get_node_record_offset(end_bound);

                    //set values for first
                    rank_in_chain=end_rank;
                    go_left_in_chain = false;
                    node_length = end_length;
                    std::tie(prefix_sum, forward_loop, reverse_loop, component) =
                            TrivialSnarlRecord (get_record_offset(end_bound), &snarl_tree_records).get_chain_values(get_node_record_offset(end_bound));

                    node_lengths_to_add += end_length;
                }
            }
        }

        if (is_node(child1) != is_node(child2) && rank_in_chain1 == rank_in_chain2 &&
               go_left1 != go_left2) {
            //If one child is a node and the other is a snarl, and the node is the boundary of the snarl
            //pointing into the snarl, then the distance is 0
            return 0;
        } else if (!is_node(child1) && !is_node(child2) && rank_in_chain1 == rank_in_chain2 &&
                go_left1 != go_left2) {
                //If the snarls are adjacent (and not the same snarl)
            return node_length2;//return the node length
        }
        

        return sum(chain_record.get_distance(rank_in_chain1, go_left1,
                                              node_length1, prefix_sum1,
                                              forward_loop1, reverse_loop1,
                                              component1, component1, 
                                              rank_in_chain2, go_left2,
                                              node_length2, prefix_sum2,
                                              forward_loop2, reverse_loop2,
                                              component2, component2),
                    node_lengths_to_add);

    } else if (is_snarl(parent)) {
        bool snarl_is_root = is_root_snarl(parent);
        size_t rank1, rank2; bool dir1, dir2;
        if (is_sentinel(child1)) {
            rank1 = starts_at(child1) == START ? 0 : 1;
            dir1 = false;
        } else {
            rank1 = get_rank_in_parent(child1);
            dir1 = !child_ends_at_start1;
        }
        if (is_sentinel(child2)) {
            rank2 = starts_at(child2) == START ? 0 : 1;
            dir2 = false;
        } else {
            rank2 = get_rank_in_parent(child2);
            dir2 = !child_ends_at_start2;
        }
        if ((is_sentinel(child1) && starts_at(child1) == ends_at(child1)) ||
            (is_sentinel(child2) && starts_at(child2) == ends_at(child2)) ) {
            //If this is a sentinel pointing out of the snarl
            return std::numeric_limits<size_t>::max();
        }


        if (get_record_type(snarl_tree_records->at(get_record_offset(parent))) == DISTANCED_SIMPLE_SNARL) {
            auto result = SimpleSnarlRecord(parent, &snarl_tree_records).get_distance(rank1, dir1, rank2, dir2);
            return result;
        } else if (is_oversized_snarl(get_record_type(snarl_tree_records->at(get_record_offset(parent))))) {
            // We need to point at the hub labeling data, which lives right after the fixed-size snarl record header
            // This points to the whole record, including the fixed-size header
            auto record_it = snarl_tree_records->begin() + get_record_offset(parent);
            // This points to the start of the packed hub label data; its first slot is the label count
            auto length_data_it = record_it + SNARL_RECORD_SIZE;
            // Here we need to convert from child rank (where 0 and 1 represent
            // snarl boundaries oriented along the snarl) and child orientation
            // to HHL vertex ranks.
            //
            // The HHL index thinks a child being "reverse" means that we're
            // thinking of the child in the opposite orientation form how it
            // appears in us. (So even not-reversed won't mean local forward
            // orientation if is_reversed_in_parent() is true for that child).
            //
            // Note that dir1 and dir2 aren't just normal is_reverse flags. 
            //
            // For a sentinel rank 1 (end node) as rank1, dir1 false needs to mean into the snarl (so start of end node, reverse strand).
            // For a sentinel rank 0 (start node) as rank1, dir1 false needs to mean into the snarl (so end of start node, forward strand).
            //
            // For a node as rank1, with its end connected to rank2, dir1 true needs to mean towards the thing attached to its end. If we're a source, that means it must be forward strand.
            //
            // For a node as rank2, with its end connected to rank1, dir2 true needs to mean towards the thing attached to its end. If we're not a source, that means it must be reverse strand.
            
            size_t from_port = bgid(rank1, !dir1 ^ (rank1 == 0), true);
            size_t to_port = bgid(rank2, dir2, false);

            size_t distance = promote_distance<size_t>(hhl_query(length_data_it, from_port, to_port));
            return distance;
          
        } else if (rank1 == 0 && rank2 == 0 && !snarl_is_root) {
            //Start to start is stored in the snarl
            auto result = SnarlRecord(parent, &snarl_tree_records).get_distance_start_start();
            return result;
        } else if ((rank1 == 0 && rank2 == 1) || (rank1 == 1 && rank2 == 0) && !snarl_is_root) {
            //start to end / end to start is stored in the snarl
            auto result = SnarlRecord(parent, &snarl_tree_records).get_min_length();
            return result;
        } else if (rank1 == 1 && rank2 == 1 && !snarl_is_root) {
            //end to end is stored in the snarl
            auto result = SnarlRecord(parent, &snarl_tree_records).get_distance_end_end();
            return result;
        } else if ((rank1 == 0 || rank1 == 1 || rank2 == 0 || rank2 == 1) && !snarl_is_root) {
            //If one node is a boundary and the other is a child
            size_t boundary_rank = (rank1 == 0 || rank1 == 1) ? rank1 : rank2;
            const net_handle_t& internal_child = (rank1 == 0 || rank1 == 1) ? child2 : child1;
            bool internal_is_reversed = (rank1 == 0 || rank1 == 1) ? dir2 : dir1;
            if (is_trivial_chain( internal_child) ) {
                //Child is just a node pretending to be a chain
                size_t result;
                if (boundary_rank == 0 && !internal_is_reversed) {
                    //Start to left of child
                    result = NodeRecord(internal_child, &snarl_tree_records).get_distance_left_start();
                } else if (boundary_rank == 0 && internal_is_reversed) {
                    //Start to right of child
                    result = NodeRecord(internal_child, &snarl_tree_records).get_distance_right_start();
                } else if (boundary_rank == 1 && !internal_is_reversed) {
                    //End to left of child
                    result = NodeRecord(internal_child, &snarl_tree_records).get_distance_left_end();
                } else {
                    //End to right of child
                    result = NodeRecord(internal_child, &snarl_tree_records).get_distance_right_end();
                }
                return result;
            } else {
                //Child is an actual chain
                size_t result;
                if (boundary_rank == 0 && !internal_is_reversed) {
                    //Start to left of child
                    result = ChainRecord(internal_child, &snarl_tree_records).get_distance_left_start();
                } else if (boundary_rank == 0 && internal_is_reversed) {
                    //Start to right of child
                    result = ChainRecord(internal_child, &snarl_tree_records).get_distance_right_start();
                } else if (boundary_rank == 1 && !internal_is_reversed) {
                    //End to left of child
                    result = ChainRecord(internal_child, &snarl_tree_records).get_distance_left_end();
                } else {
                    //End to right of child
                    result = ChainRecord(internal_child, &snarl_tree_records).get_distance_right_end();
                }
                return result;
            }
        } else {
            auto result = SnarlRecord(parent, &snarl_tree_records).get_distance(rank1, dir1, rank2, dir2);
            return result;
        }
    } else {
        throw runtime_error("error: Trying to find distance in the wrong type of handle");
    }
}

size_t SnarlDistanceIndex::distance_in_snarl(const net_handle_t& parent, 
        const size_t& rank1, const bool& right_side1, const size_t& rank2, const bool& right_side2, 
        const HandleGraph* graph, size_t distance_limit) const {

    bool snarl_is_root = is_root_snarl(parent);
    
    if (get_record_type(snarl_tree_records->at(get_record_offset(parent))) == DISTANCED_SIMPLE_SNARL) {
        return SimpleSnarlRecord(parent, &snarl_tree_records).get_distance(rank1, right_side1, rank2, right_side2);
    } else if (is_oversized_snarl(get_record_type(snarl_tree_records->at(get_record_offset(parent))))
        && !(rank1 == 0 || rank1 == 1 || rank2 == 0 || rank2 == 1) ) {
        //If this is an oversized snarl and we're looking for internal distances, use the hub labels.
        auto record_it = snarl_tree_records->begin() + get_record_offset(parent);
        auto length_data_it = record_it + SNARL_RECORD_SIZE;
        size_t from_port = bgid(rank1, !right_side1 ^ (rank1 == 0), true);
        size_t to_port   = bgid(rank2,  right_side2,                false);
        return promote_distance<size_t>(hhl_query(length_data_it, from_port, to_port));

    } else if (rank1 == 0 && rank2 == 0 && !snarl_is_root) {
        //Start to start is stored in the snarl
        return SnarlRecord(parent, &snarl_tree_records).get_distance_start_start();
    } else if (((rank1 == 0 && rank2 == 1) || (rank1 == 1 && rank2 == 0)) && !snarl_is_root) {
        //start to end / end to start is stored in the snarl
        return SnarlRecord(parent, &snarl_tree_records).get_min_length();
    } else if (rank1 == 1 && rank2 == 1 && !snarl_is_root) {
        //end to end is stored in the snarl
        return SnarlRecord(parent, &snarl_tree_records).get_distance_end_end();
    } else if ((rank1 == 0 || rank1 == 1 || rank2 == 0 || rank2 == 1) && !snarl_is_root) {
        //If one node is a boundary and the other is a child
        size_t boundary_rank = (rank1 == 0 || rank1 == 1) ? rank1 : rank2;
        const net_handle_t& internal_child = (rank1 == 0 || rank1 == 1) ? get_snarl_child_from_rank(parent, rank2) 
                                                                        : get_snarl_child_from_rank(parent, rank1);
        bool internal_is_reversed = (rank1 == 0 || rank1 == 1) ? right_side2 : right_side1;
        if (is_trivial_chain( internal_child) ) {
            //Child is just a node pretending to be a chain
            if (boundary_rank == 0 && !internal_is_reversed) {
                //Start to left of child
                return NodeRecord(internal_child, &snarl_tree_records).get_distance_left_start();
            } else if (boundary_rank == 0 && internal_is_reversed) {
                //Start to right of child
                return NodeRecord(internal_child, &snarl_tree_records).get_distance_right_start();
            } else if (boundary_rank == 1 && !internal_is_reversed) {
                //End to left of child
                return NodeRecord(internal_child, &snarl_tree_records).get_distance_left_end();
            } else {
                //End to right of child
                return NodeRecord(internal_child, &snarl_tree_records).get_distance_right_end();
            }
        } else {
            //Child is an actual chain
            if (boundary_rank == 0 && !internal_is_reversed) {
                //Start to left of child
                return ChainRecord(internal_child, &snarl_tree_records).get_distance_left_start();
            } else if (boundary_rank == 0 && internal_is_reversed) {
                //Start to right of child
                return ChainRecord(internal_child, &snarl_tree_records).get_distance_right_start();
            } else if (boundary_rank == 1 && !internal_is_reversed) {
                //End to left of child
                return ChainRecord(internal_child, &snarl_tree_records).get_distance_left_end();
            } else {
                //End to right of child
                return ChainRecord(internal_child, &snarl_tree_records).get_distance_right_end();
            }
        }
    } else {
       return SnarlRecord(get_record_offset(parent), &snarl_tree_records).get_distance(rank1, right_side1, rank2, right_side2);
    }

}
size_t SnarlDistanceIndex::max_distance_in_parent(const net_handle_t& parent, 
        const net_handle_t& child1, const net_handle_t& child2, const HandleGraph* graph, size_t distance_limit) const {


    if (is_chain(parent)) {

        //Get the orientation of the children. This only cares about the end endpoint, and assumes that things that end
        //in tips are start-end, meaning that they don't end at the start
        bool child_ends_at_start1 = ends_at(child1) == START;
        bool child_ends_at_start2 = ends_at(child2) == START;

        if (get_record_handle_type(get_record_type(snarl_tree_records->at(get_record_offset(parent)))) == NODE_HANDLE ||
            get_record_handle_type(get_record_type(snarl_tree_records->at(get_record_offset(parent)))) == SNARL_HANDLE) {
            return std::numeric_limits<size_t>::max();
        }
        ChainRecord chain_record(parent, &snarl_tree_records);

        //These are the values needed for calculating the distance in the chain
        //They are:
        //rank in chain, go_left, node_length, prefix_sum, forward_loop, reverse_loop, component
        size_t rank_in_chain1; size_t rank_in_chain2;
        bool go_left1; bool go_left2;
        size_t node_length1; size_t node_length2;
        size_t prefix_sum1; size_t prefix_sum2;
        size_t forward_loop1 ; size_t forward_loop2;
        size_t reverse_loop1 ; size_t reverse_loop2;
        size_t component1; size_t component2;
        size_t end_component1; size_t end_component2;

        //The node lengths of snarl boundaries, to be added to the total distance
        size_t node_lengths_to_add = 0;

        //Get all the relevant chain values
        for (size_t i = 0 ; i <=1 ; i ++ ) {
            const net_handle_t& child = i == 0 ? child1 : child2;
            size_t& rank_in_chain= i == 0 ? rank_in_chain1: rank_in_chain2;
            bool& go_left_in_chain = i == 0 ? go_left1  : go_left2;
            size_t& node_length  = i == 0 ? node_length1  : node_length2;
            size_t& prefix_sum   = i == 0 ? prefix_sum1   : prefix_sum2;
            size_t& forward_loop = i == 0 ? forward_loop1 : forward_loop2;
            size_t& reverse_loop = i == 0 ? reverse_loop1 : reverse_loop2;
            size_t& component    = i == 0 ? component1    : component2;
            size_t& end_component    = i == 0 ? end_component1    : end_component2;
            bool& child_ends_at_start = i == 0 ? child_ends_at_start1 : child_ends_at_start2;

            if (is_node(child)) {

                bool go_left = is_reversed_in_parent(child);
                if (child_ends_at_start){
                    go_left = !go_left;
                }
                rank_in_chain=get_rank_in_parent(child);
                go_left_in_chain = go_left;
                node_length = minimum_length(child);
                std::tie(prefix_sum, forward_loop, reverse_loop, component) = 
                    TrivialSnarlRecord(get_record_offset(child), &snarl_tree_records).get_chain_values(get_node_record_offset(child));
                prefix_sum = TrivialSnarlRecord(get_record_offset(child), &snarl_tree_records).get_max_prefix_sum(get_node_record_offset(child));
                end_component = TrivialSnarlRecord(get_record_offset(child), &snarl_tree_records).get_chain_component(get_node_record_offset(child), true);

                
            } else {
                //If this is a snarl, then find the appropriate boundary node and fill it in
                assert(is_snarl(child));
                if (child_ends_at_start) {
                    //If we're going backwards in the chain, then we want the start node
                    //And turn it into a node, since it will be a sentinel
                    net_handle_t start_bound = get_node_from_sentinel(get_bound(child, false, false));

                    size_t start_length = minimum_length(start_bound);

                        //set values for first
                    rank_in_chain=get_record_offset(start_bound) + get_node_record_offset(start_bound);
                    go_left_in_chain = true;
                    node_length = start_length;
                    std::tie(prefix_sum, forward_loop, reverse_loop, component) = 
                        TrivialSnarlRecord (get_record_offset(start_bound), &snarl_tree_records).get_chain_values(get_node_record_offset(start_bound));

                    prefix_sum = TrivialSnarlRecord (get_record_offset(start_bound), &snarl_tree_records).get_max_prefix_sum(get_node_record_offset(start_bound));
                    node_lengths_to_add += start_length;
                } else {
                    //Do the same thing for the snarl end node if we're going forwards
                    net_handle_t end_bound = get_node_from_sentinel(get_bound(child, true, false));
                    size_t end_length = minimum_length(end_bound);
                    size_t end_rank = get_record_offset(end_bound) + get_node_record_offset(end_bound);

                    //set values for first
                    rank_in_chain=end_rank;
                    go_left_in_chain = false;
                    node_length = end_length;
                    std::tie(prefix_sum, forward_loop, reverse_loop, component) =
                            TrivialSnarlRecord (get_record_offset(end_bound), &snarl_tree_records).get_chain_values(get_node_record_offset(end_bound));

                    prefix_sum = TrivialSnarlRecord (get_record_offset(end_bound), &snarl_tree_records).get_max_prefix_sum(get_node_record_offset(end_bound));

                    node_lengths_to_add += end_length;
                }
            }
        }

        if (is_node(child1) != is_node(child2) && rank_in_chain1 == rank_in_chain2 &&
               go_left1 != go_left2) {
            //If one child is a node and the other is a snarl, and the node is the boundary of the snarl
            //pointing into the snarl, then the distance is 0
            return 0;
        } else if (!is_node(child1) && !is_node(child2) && rank_in_chain1 == rank_in_chain2 &&
                go_left1 != go_left2) {
                //If the snarls are adjacent (and not the same snarl)
            return node_length2;//return the node length
        }

        return sum(chain_record.get_distance(rank_in_chain1, go_left1,
                                              node_length1, prefix_sum1,
                                              forward_loop1, reverse_loop1,
                                              component1, component1, 
                                              rank_in_chain2, go_left2,
                                              node_length2, prefix_sum2,
                                              forward_loop2, reverse_loop2,
                                              component2, component2),
                    node_lengths_to_add);

    } else {
        //If the parent isn't a chain, then just return the minimum
        return distance_in_parent(parent, child1, child2, graph, distance_limit);
    }
}

size_t SnarlDistanceIndex::distance_to_parent_bound(const net_handle_t& parent, bool to_start, net_handle_t child,
    tuple<net_handle_record_t, net_handle_record_t, net_handle_record_t, net_handle_record_t> parent_and_child_types) const {
    /* parent_and_child_types is a tuple of parent handle type, parent record type, child handle type, child record type
      * This is really just used to see if the parent and child are trivial chains, so it might not be exactly what 
      * the actual record is
      * */



    //Get the record and handle types of the parent and child
    bool has_handle_types = (parent_and_child_types != make_tuple(ROOT_HANDLE, ROOT_HANDLE, ROOT_HANDLE, ROOT_HANDLE));
    bool parent_is_chain = has_handle_types ? std::get<0>(parent_and_child_types) == CHAIN_HANDLE 
                                            : is_chain(parent);
    bool parent_is_trivial_chain = has_handle_types ? (std::get<0>(parent_and_child_types) == NODE_HANDLE && std::get<1>(parent_and_child_types) == CHAIN_HANDLE)
                                                    : is_trivial_chain(parent);
    bool child_is_trivial_chain = has_handle_types ? (std::get<2>(parent_and_child_types) == NODE_HANDLE && std::get<3>(parent_and_child_types) == CHAIN_HANDLE)
                                                    : is_trivial_chain(child);
    bool parent_is_snarl = has_handle_types ? std::get<0>(parent_and_child_types) == SNARL_HANDLE
                                            : is_snarl(parent);

    //Get the orientation of the child. This only cares about the end endpoint, and assumes that things that end
    //in tips are start-end, meaning that they don't end at the start
    bool child_ends_at_start = ends_at(child) == START;
    if (parent_is_trivial_chain) {
        //If this is a node pretending to be a chain and the parent is the trivial chain, 
        //then the child is always oriented forwards in the parent (because it is the parent) and the distance is 0 or inf
        if ((child_ends_at_start && to_start) ||
            (!child_ends_at_start && !to_start)) {
            //If the child is traversed backwards and we're going to the start
            //Or the child is traversed forwards and we're going to the end
            return 0;
        } else {
            return std::numeric_limits<size_t>::max();
        } 
    } else if (parent_is_snarl) {
        if (is_simple_snarl(parent)) {
            if (to_start) {
                if (child_ends_at_start != is_reversed_in_parent(child)) {
                    //If the child is going backwards, and it isn't reversed in the parent,
                    //or if the child is going forwards and it is reversed in the parent,
                    //then the distances to the start will be 0
                    return 0;
                } else {
                    return std::numeric_limits<size_t>::max();
                }
            } else {
                //to end of snarl
                if (child_ends_at_start == is_reversed_in_parent(child)) {
                    //If the child is going backwards and it is reversed in the parent
                    //or the child is going forwards and it isn't reversed in the parent
                    return 0;
                } else {
                    return std::numeric_limits<size_t>::max();
                }
            }
        } else {
            //Otherwise, the parent is a normal snarl and the child nodes/chains know the distance
            //to the ends of the snarl
    
            if (child_is_trivial_chain) {
                NodeRecord child_record (child, &snarl_tree_records);
                if (child_ends_at_start && to_start) {
                    return child_record.get_distance_left_start();
                } else if (!child_ends_at_start && to_start) {
                    return child_record.get_distance_right_start();
                } else if (child_ends_at_start && !to_start) { 
                    return child_record.get_distance_left_end();
                } else {
                    return child_record.get_distance_right_end();
                }
            } else {
                ChainRecord child_record (child, &snarl_tree_records);
                if (child_ends_at_start && to_start) {
                    return child_record.get_distance_left_start();
                } else if (!child_ends_at_start && to_start) {
                    return child_record.get_distance_right_start();
                } else if (child_ends_at_start && !to_start) { 
                    return child_record.get_distance_left_end();
                } else {
                    return child_record.get_distance_right_end();
                }
            }
        }


    } else {
        //Otherwise, the parent is a chain so get the bound and use distance_in_parent


        //The parent bound we want, facing in
        net_handle_t parent_bound = get_bound(parent, !to_start, true);

        if (parent_bound == flip(child)) {
            //If the child is the bound we want and pointing into it (out of the chain)
            return 0;
        }
        //The node length of the boundary node, only set for chains 
        size_t bound_length = !parent_is_chain ? 0 : minimum_length(parent_bound);

        if (parent_is_chain && !to_start){
            if (is_looping_chain(parent)){
                //If this is a looping chain and we want the end 
                bound_length = 0;
            }
        }
        
        return sum(bound_length, 
                    distance_in_parent(parent, parent_bound, child));
    }

}

tuple<nid_t, bool, bool> SnarlDistanceIndex::into_which_snarl(const nid_t& id, const bool& reverse) const {
    net_handle_t node = get_node_net_handle(id, reverse);
    net_handle_t parent = get_parent(node);

    if (is_trivial_chain(parent) || is_root(parent)) {
        //If the parent is a trivial chain, then it is a node child of a snarl
        return {0, false, false};
    } else {
        //Otherwise, the parent is an actual chain and this is the boundary node of a (possibly trivial) snarl

        //If the node is traversed forwards in the parent and we are going backwards (or the opposite),
        //then we are going backwards in the chain
        bool go_left_in_chain = is_reversed_in_parent(node) != reverse;
        ChainRecord chain_record (parent, &snarl_tree_records);
        net_handle_t next_child = chain_record.get_next_child(node, go_left_in_chain);
        if (is_snarl(next_child)) {
            //If the node points into a snarl
            net_handle_t start = get_node_from_sentinel(get_bound(next_child, false, true));
            return {node_id(start), ends_at(start) == START, false}; 
        } else if (next_child == node) {
            //If the node points out of the chain
            return {0, false, false};
        } else {
            //If the node points to another node (trivial snarl)
            //Assume that the trivial snarl is oriented forwards in the chain, so the start node
            //is whichever comes first in the chain
            if (go_left_in_chain) {
                //If we went left in the chain, then the start node of the snarl is the next node
                return {node_id(next_child), is_reversed_in_parent(next_child), true};
            } else {
                //If we went right in the chain, then the original node is the start node
                return {id, reverse, true};
            }
        }
    }
}




pair<net_handle_t, bool> SnarlDistanceIndex::lowest_common_ancestor(const net_handle_t& net1, const net_handle_t& net2) const {
    net_handle_t parent1 = net1;
    net_handle_t parent2 = net2;

    std::unordered_set<net_handle_t> net1_ancestors;
    while (!is_root(parent1)){
        net1_ancestors.insert(canonical(parent1));
        parent1 = get_parent(parent1);
    }

    while (net1_ancestors.count(canonical(parent2)) == 0 && !is_root(parent2)){
        //Go up until the parent2 matches something in the ancestors of net1
        //This loop will end because everything is in the same root eventually
        parent2 = get_parent(parent2);
    }

    bool is_connected = true;
    if (is_root(parent2) && is_root(parent1)){
        size_t parent_record_offset1 = SnarlTreeRecord(parent1, &snarl_tree_records).get_parent_record_offset();
        size_t parent_record_offset2 = SnarlTreeRecord(parent2, &snarl_tree_records).get_parent_record_offset();
        if (parent_record_offset1 != parent_record_offset2) {
            is_connected = false;
        }
    }
    return make_pair(canonical(parent2), is_connected);
}

size_t SnarlDistanceIndex::minimum_distance(const handlegraph::nid_t id1, const bool rev1, const size_t offset1, 
                                            const handlegraph::nid_t id2, const bool rev2, const size_t offset2, 
                                            bool unoriented_distance, const HandleGraph* graph) const {
    RootRecord root_record (get_root(), &snarl_tree_records);
    size_t max_node_id = root_record.get_min_node_id() + root_record.get_node_count();
    if (id1 < root_record.get_min_node_id() || id2 < root_record.get_min_node_id() ||
        id1 > max_node_id || id2 > max_node_id) {
        throw runtime_error("error: Looking for the minimum distance of a node that does not exist");
    }




    /*Helper function to walk up the snarl tree
     * Given a net handle, its parent,  and the distances to the start and end of the handle, 
     * update the distances to reach the ends of the parent and update the handle and its parent
     * If the parent is a chain, then the new distances include the boundary nodes of the chain.
     * If it is a snarl, it does not*/
    auto update_distances = [&](net_handle_t& net, net_handle_t& parent, size_t& dist_start, size_t& dist_end, bool first_node) {

        if (is_trivial_chain(parent)) {
            //Don't update distances for the trivial chain
            return;
        } else if (is_simple_snarl(parent)) {
            //If it's a simple snarl just check if they should be reversed
            if (is_reversed_in_parent (net)) {
                size_t tmp = dist_start;
                dist_start = dist_end;
                dist_end = tmp;
            } 

            return;
        }


        net_handle_t start_bound = get_bound(parent, false, true);
        net_handle_t end_bound = get_bound(parent, true, true);

        //The lengths of the start and end nodes of net
        //This is only needed if net is a snarl, since the boundary nodes are not technically part of the snarl
        size_t start_length = is_chain(parent) ? node_length(start_bound) : 0;
        size_t end_length = is_chain(parent) ? node_length(end_bound) : 0;

        //Get the distances from the bounds of the parent to the node we're looking at
        size_t distance_start_start = start_bound == net ? 0 
                : sum(start_length, distance_in_parent(parent, start_bound, flip(net), graph));
        size_t distance_start_end = start_bound == flip(net) ? 0 
                : sum(start_length, distance_in_parent(parent, start_bound, net, graph));
        size_t distance_end_start = end_bound == net ? 0 
                : sum(end_length, distance_in_parent(parent, end_bound, flip(net), graph));
        size_t distance_end_end = end_bound == flip(net) ? 0 
                : sum(end_length, distance_in_parent(parent, end_bound, net, graph));

        size_t distance_start = dist_start;
        size_t distance_end = dist_end; 



        dist_start = std::min( sum(distance_start_start, distance_start), 
                                    sum(distance_start_end , distance_end));
        dist_end = std::min(sum(distance_end_start , distance_start), 
                            sum(distance_end_end , distance_end));
        return;
    };

    /*
     * Get net handles for the two nodes and the distances from each position to the ends of the handles
     */
    net_handle_t net1 = get_node_net_handle(id1);
    net_handle_t net2 = get_node_net_handle(id2);
    pair<net_handle_t, bool> lowest_ancestor = lowest_common_ancestor(net1, net2);
    if (!lowest_ancestor.second) {
        //If these are not in the same connected component
        return std::numeric_limits<size_t>::max();
    }

    //The lowest common ancestor of the two positions
    net_handle_t common_ancestor = start_end_traversal_of(lowest_ancestor.first);

    //These are the distances to the ends of the node, including the position
    size_t distance_to_start1 = rev1 ? node_length(net1) - offset1 : offset1 + 1;
    size_t distance_to_end1 = rev1 ? offset1 + 1 : node_length(net1) - offset1;
    size_t distance_to_start2 = rev2 ? node_length(net2) - offset2 : offset2 + 1;
    size_t distance_to_end2 = rev2 ? offset2 + 1 : node_length(net2) - offset2;

    if (!unoriented_distance) {
        //If we care about the oriented distance, one of the distances will be infinite
        if (rev1) {
            distance_to_end1 = std::numeric_limits<size_t>::max();
        } else {
            distance_to_start1 = std::numeric_limits<size_t>::max();
        }
        if (rev2) {
            distance_to_start2 = std::numeric_limits<size_t>::max();
        } else {
            distance_to_end2 = std::numeric_limits<size_t>::max();
        }
    }


    size_t minimum_distance = std::numeric_limits<size_t>::max();



    /*
     * Walk up the snarl tree until net1 and net2 are children of the lowest common ancestor
     * Keep track of the distances to the ends of the net handles as we go
     */
 
    if (start_end_traversal_of(net1) == start_end_traversal_of(net2)){
        if (sum(distance_to_end1 , distance_to_start2) > node_length(net1) && 
            sum(distance_to_end1 , distance_to_start2) != std::numeric_limits<size_t>::max()) {
            //If the positions are on the same node and are pointing towards each other, then
            //check the distance between them in the node
            minimum_distance = minus(sum(distance_to_end1 , distance_to_start2), node_length(net1));
        }
        if (sum(distance_to_start1 , distance_to_end2) > node_length(net1) && 
            sum(distance_to_start1 , distance_to_end2) != std::numeric_limits<size_t>::max()) {
            minimum_distance = std::min(minus(sum(distance_to_start1 , distance_to_end2), node_length(net1)), minimum_distance);
        }
        common_ancestor = start_end_traversal_of(get_parent(net1));
    } else {

        //Get the distance from position 1 up to the ends of a child of the common ancestor
        while (start_end_traversal_of(get_parent(net1)) != common_ancestor && !is_root(get_parent(net1))) {
            net_handle_t parent = start_end_traversal_of(get_parent(net1));
            update_distances(net1, parent, distance_to_start1, distance_to_end1, true);
            net1 = parent;
        }
        //And the same for position 2
        while (start_end_traversal_of(get_parent(net2)) != start_end_traversal_of(common_ancestor) && !is_root(get_parent(net2))) {
            net_handle_t parent = start_end_traversal_of(get_parent(net2));
            update_distances(net2, parent, distance_to_start2, distance_to_end2, false);
            net2 = parent;
        }
    }


    /* 
     * common_ancestor is now the lowest common ancestor of both net handles, and 
     * net1 and net2 are both children of common_ancestor
     * Walk up to the root and check for distances between the positions within each
     * ancestor
     */

    while (!is_root(net1)){

        //Find the minimum distance between the two children (net1 and net2).
        //Since distance_in_parent() takes traversals that point toward each
        //other, and not in the same direction along the connecting path,
        //distance_end_start (the "normal" forward/forward case) has to flip
        //net2 to point back towards net1, and the others behave similarly.
        size_t distance_start_start = distance_in_parent(common_ancestor, flip(net1), flip(net2), graph);
        size_t distance_start_end = distance_in_parent(common_ancestor, flip(net1), net2, graph);
        size_t distance_end_start = distance_in_parent(common_ancestor, net1, flip(net2), graph);
        size_t distance_end_end = distance_in_parent(common_ancestor, net1, net2, graph);


        //And add those to the distances we've found to get the minimum distance between the positions
        minimum_distance = std::min(minimum_distance, 
                           std::min(sum(sum(distance_start_start , distance_to_start1), distance_to_start2),
                           std::min(sum(sum(distance_start_end , distance_to_start1), distance_to_end2),
                           std::min(sum(sum(distance_end_start , distance_to_end1), distance_to_start2),
                                    sum(sum(distance_end_end , distance_to_end1), distance_to_end2)))));


        if (!is_root(common_ancestor)) {
            //Update the distances to reach the ends of the common ancestor
            update_distances(net1, common_ancestor, distance_to_start1, distance_to_end1, true);
            update_distances(net2, common_ancestor, distance_to_start2, distance_to_end2, false);

            //Update which net handles we're looking at
            net1 = common_ancestor;
            net2 = common_ancestor;
            common_ancestor = start_end_traversal_of(get_parent(common_ancestor));
        } else {
            //Just update this one to break out of the loop
            net1 = common_ancestor;

        }

    }


    //minimum distance currently includes both positions
    return minimum_distance == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() : minimum_distance-1;



}
size_t SnarlDistanceIndex::maximum_distance(const handlegraph::nid_t id1, const bool rev1, const size_t offset1, 
                                            const handlegraph::nid_t id2, const bool rev2, const size_t offset2, 
                                            bool unoriented_distance, const HandleGraph* graph) const {
    RootRecord root_record (get_root(), &snarl_tree_records);
    size_t max_node_id = root_record.get_min_node_id() + root_record.get_node_count();
    if (id1 < root_record.get_min_node_id() || id2 < root_record.get_min_node_id() ||
        id1 > max_node_id || id2 > max_node_id) {
        throw runtime_error("error: Looking for the maximum distance of a node that does not exist");
    }




    /*Helper function to walk up the snarl tree
     * Given a net handle, its parent,  and the distances to the start and end of the handle, 
     * update the distances to reach the ends of the parent and update the handle and its parent
     * If the parent is a chain, then the new distances include the boundary nodes of the chain.
     * If it is a snarl, it does not*/
    //TODO: This should really be an actual function and it doesn't consider the lengths of the 
    //boundary nodes. I think chains might need to know the lengths of their boundary nodes,
    //or it could find it from the snarl. Really, since the snarl is storing it anyway, I should
    //probaby rearrange things so that either the snarl or the chain can access boundary node lengths
    auto update_distances = [&](net_handle_t& net, net_handle_t& parent, size_t& dist_start, size_t& dist_end, bool first_node) {

        if (is_trivial_chain(parent)) {
            //Don't update distances for the trivial chain
            return;
        } else if (is_simple_snarl(parent)) {
            //If it's a simple snarl just check if they should be reversed
            if (is_reversed_in_parent (net)) {
                size_t tmp = dist_start;
                dist_start = dist_end;
                dist_end = tmp;
            } 

            return;
        }


        net_handle_t start_bound = get_bound(parent, false, true);
        net_handle_t end_bound = get_bound(parent, true, true);

        //The lengths of the start and end nodes of net
        //This is only needed if net is a snarl, since the boundary nodes are not technically part of the snarl
        size_t start_length = is_chain(parent) ? node_length(start_bound) : 0;
        size_t end_length = is_chain(parent) ? node_length(end_bound) : 0;

        //Get the distances from the bounds of the parent to the node we're looking at
        size_t distance_start_start = start_bound == net ? 0 
                : sum(start_length, max_distance_in_parent(parent, start_bound, flip(net), graph));
        size_t distance_start_end = start_bound == flip(net) ? 0 
                : sum(start_length, max_distance_in_parent(parent, start_bound, net, graph));
        size_t distance_end_start = end_bound == net ? 0 
                : sum(end_length, max_distance_in_parent(parent, end_bound, flip(net), graph));
        size_t distance_end_end = end_bound == flip(net) ? 0 
                : sum(end_length, max_distance_in_parent(parent, end_bound, net, graph));

        size_t distance_start = dist_start;
        size_t distance_end = dist_end; 

        dist_start = maximum( sum(distance_start_start, distance_start), 
                              sum(distance_start_end , distance_end));
        dist_end = maximum(sum(distance_end_start , distance_start), 
                           sum(distance_end_end , distance_end));
        return;
    };

    /*
     * Get net handles for the two nodes and the distances from each position to the ends of the handles
     * TODO: net2 is pointing in the opposite direction of the position. The final 
     * distance will be between the two nets pointing towards each other
     */
    net_handle_t net1 = get_node_net_handle(id1);
    net_handle_t net2 = get_node_net_handle(id2);
    pair<net_handle_t, bool> lowest_ancestor = lowest_common_ancestor(net1, net2);
    if (!lowest_ancestor.second) {
        //If these are not in the same connected component
        return std::numeric_limits<size_t>::max();
    }

    //The lowest common ancestor of the two positions
    net_handle_t common_ancestor = start_end_traversal_of(lowest_ancestor.first);

    //These are the distances to the ends of the node, including the position
    size_t distance_to_start1 = rev1 ? node_length(net1) - offset1 : offset1 + 1;
    size_t distance_to_end1 = rev1 ? offset1 + 1 : node_length(net1) - offset1;
    size_t distance_to_start2 = rev2 ? node_length(net2) - offset2 : offset2 + 1;
    size_t distance_to_end2 = rev2 ? offset2 + 1 : node_length(net2) - offset2;

    if (!unoriented_distance) {
        //If we care about the oriented distance, one of the distances will be infinite
        if (rev1) {
            distance_to_end1 = std::numeric_limits<size_t>::max();
        } else {
            distance_to_start1 = std::numeric_limits<size_t>::max();
        }
        if (rev2) {
            distance_to_start2 = std::numeric_limits<size_t>::max();
        } else {
            distance_to_end2 = std::numeric_limits<size_t>::max();
        }
    }


    size_t maximum_distance = std::numeric_limits<size_t>::max();

    //Remember which common ancestor the maximum distance actually comes from
    //tuple of <the ancestor, distance between the children, connectivity between the two children
    tuple<net_handle_t, size_t, connectivity_t> common_ancestor_connectivity (get_root(), std::numeric_limits<size_t>::max(), START_END); 


    /*
     * Walk up the snarl tree until net1 and net2 are children of the lowest common ancestor
     * Keep track of the distances to the ends of the net handles as we go
     */
 
    if (start_end_traversal_of(net1) == start_end_traversal_of(net2)){
        if (sum(distance_to_end1 , distance_to_start2) > node_length(net1) && 
            sum(distance_to_end1 , distance_to_start2) != std::numeric_limits<size_t>::max()) {
            //If the positions are on the same node and are pointing towards each other, then
            //check the distance between them in the node
            maximum_distance = minus(sum(distance_to_end1 , distance_to_start2), node_length(net1));
        }
        if (sum(distance_to_start1 , distance_to_end2) > node_length(net1) && 
            sum(distance_to_start1 , distance_to_end2) != std::numeric_limits<size_t>::max()) {
            size_t new_max = minus(sum(distance_to_start1 , distance_to_end2), node_length(net1));
            if (new_max != std::numeric_limits<size_t>::max()) {
                maximum_distance = maximum_distance == std::numeric_limits<size_t>::max() ? new_max 
                                 : std::max(new_max, maximum_distance);
            }
        }
        common_ancestor = start_end_traversal_of(get_parent(net1));
    } else {

        //Get the distance from position 1 up to the ends of a child of the common ancestor
        while (start_end_traversal_of(get_parent(net1)) != common_ancestor && !is_root(get_parent(net1))) {
            net_handle_t parent = start_end_traversal_of(get_parent(net1));
            update_distances(net1, parent, distance_to_start1, distance_to_end1, true);
            net1 = parent;
        }
        //And the same for position 2
        while (start_end_traversal_of(get_parent(net2)) != start_end_traversal_of(common_ancestor) && !is_root(get_parent(net2))) {
            net_handle_t parent = start_end_traversal_of(get_parent(net2));
            update_distances(net2, parent, distance_to_start2, distance_to_end2, false);
            net2 = parent;
        }
    }
    //TODO: I'm taking this out because it should really be start-end connected, but this
    //won't do this if that traversal isn't possible. Really it should just be setting the 
    //connectivity instead
    //net1 = canonical(net1);
    //net2 = canonical(net2);


    /* 
     * common_ancestor is now the lowest common ancestor of both net handles, and 
     * net1 and net2 are both children of common_ancestor
     * Walk up to the root and check for distances between the positions within each
     * ancestor
     */

    while (!is_root(net1)){
        //TODO: Actually checking distance in a chain is between the nodes, not the snarl
        //and neither include the lengths of the nodes

        //Find the maximum distance between the two children (net1 and net2)
        size_t distance_start_start = max_distance_in_parent(common_ancestor, flip(net1), flip(net2), graph);
        size_t distance_start_end = max_distance_in_parent(common_ancestor, flip(net1), net2, graph);
        size_t distance_end_start = max_distance_in_parent(common_ancestor, net1, flip(net2), graph);
        size_t distance_end_end = max_distance_in_parent(common_ancestor, net1, net2, graph);

        size_t old_maximum = maximum_distance;

        //And add those to the distances we've found to get the maximum distance between the positions
        maximum_distance = maximum(maximum_distance, 
                           maximum(sum(sum(distance_start_start , distance_to_start1), distance_to_start2),
                           maximum(sum(sum(distance_start_end , distance_to_start1), distance_to_end2),
                           maximum(sum(sum(distance_end_start , distance_to_end1), distance_to_start2),
                                    sum(sum(distance_end_end , distance_to_end1), distance_to_end2)))));


        if (!is_root(common_ancestor)) {
            //Update the distances to reach the ends of the common ancestor
            update_distances(net1, common_ancestor, distance_to_start1, distance_to_end1, true);
            update_distances(net2, common_ancestor, distance_to_start2, distance_to_end2, false);

            //Update which net handles we're looking at
            net1 = common_ancestor;
            net2 = common_ancestor;
            common_ancestor = start_end_traversal_of(get_parent(common_ancestor));
        } else {
            //Just update this one to break out of the loop
            net1 = common_ancestor;

        }

    }


    //maximum distance currently includes both positions
    return maximum_distance == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() : maximum_distance-1;

}


size_t SnarlDistanceIndex::node_length(const net_handle_t& net) const {
    if (is_node(net)) {
        auto record_type = get_record_type(snarl_tree_records->at(get_record_offset(net)));
        if (record_type == DISTANCED_NODE) {
            return NodeRecord(net, &snarl_tree_records).get_node_length();
        } else if (record_type == DISTANCED_SIMPLE_SNARL) {
            return SimpleSnarlRecord(net, &snarl_tree_records).get_node_length();
        } else {
            assert(record_type == DISTANCED_TRIVIAL_SNARL);
            return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_node_length(get_node_record_offset(net));
        }
    } else if (is_sentinel(net)) {
        return node_length(get_node_from_sentinel(net));
    } else {
        throw runtime_error("error: Looking for the node length of a non-node net_handle_t");
    }

}


size_t SnarlDistanceIndex::minimum_length(const net_handle_t& net) const {
    auto record_type = SnarlTreeRecord(net, &snarl_tree_records).get_record_type();
    if (!has_distances(record_type)) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    } else if (record_type == DISTANCED_TRIVIAL_SNARL) {
        return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_node_length(get_node_record_offset(net));
    } else if (record_type == DISTANCED_SIMPLE_SNARL) {
        if (is_snarl(net)) {
            return SimpleSnarlRecord(net, &snarl_tree_records).get_min_length();
        } else {
            return SimpleSnarlRecord(net, &snarl_tree_records).get_node_length();
        }
    } else if (record_type == MULTICOMPONENT_CHAIN) {
        return std::numeric_limits<size_t>::max();
    } else {
        return SnarlTreeRecord(net, &snarl_tree_records).get_min_length();
    }
}

size_t SnarlDistanceIndex::chain_minimum_length(const net_handle_t& net) const {
    if (is_trivial_chain(net)) {
        return minimum_length(net);
    } else {
        return ChainRecord(net, &snarl_tree_records).get_min_length();
    }
}
size_t SnarlDistanceIndex::maximum_length(const net_handle_t& net) const {
    auto record_type = SnarlTreeRecord(net, &snarl_tree_records).get_record_type();
    if (!has_distances(record_type)) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    } else if (record_type == DISTANCED_TRIVIAL_SNARL &&
        get_node_record_offset(net) != 0) {
        return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_node_length(get_node_record_offset(net));
    } else if (record_type == DISTANCED_SIMPLE_SNARL) {
        if (is_snarl(net)) {
            return SimpleSnarlRecord(net, &snarl_tree_records).get_max_length();
        } else {
            return SimpleSnarlRecord(net, &snarl_tree_records).get_node_length();
        }
    } else {
        return SnarlTreeRecord(net, &snarl_tree_records).get_max_length();
    }
}
nid_t SnarlDistanceIndex::node_id(const net_handle_t& net) const {
    if (is_node(net) || is_trivial_chain(net)) {
        if (is_node(get_record_type(snarl_tree_records->at(get_record_offset(net))))) {
            return NodeRecord(net, &snarl_tree_records).get_node_id();
        }  else if (is_simple_snarl(get_record_type(snarl_tree_records->at(get_record_offset(net))))) {
            return SimpleSnarlRecord(net, &snarl_tree_records).get_node_id();
        } else {
            return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_node_id(get_node_record_offset(net));
        }
    } else if (is_sentinel(net)) {
        // We don't know if this is a trivial or nontrivial snarl, so we need
        // to access it with the base class.
        // TODO: Make a base class for any kind of snarl.
        SnarlTreeRecord snarl_record(net, &snarl_tree_records);
        if (get_start_endpoint(net) == START) {
            return snarl_record.get_start_id();
        } else {
            return snarl_record.get_end_id();
        }
    } else {
        throw runtime_error("error: Looking for the node id of a non-node net_handle_t");
    }

}
bool SnarlDistanceIndex::has_node(const nid_t id) const {
    RootRecord root_record (get_root(), &snarl_tree_records);
    size_t node_pointer_offset = get_node_pointer_offset(id, root_record.get_min_node_id(), root_record.get_connected_component_count());
    return snarl_tree_records->at(node_pointer_offset) != 0;
}

bool SnarlDistanceIndex::is_reversed_in_parent(const net_handle_t& net) const {
    SnarlTreeRecord record(net, &snarl_tree_records);
    if (is_trivial_snarl(record.get_record_type())) {
        return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_is_reversed_in_parent(get_node_record_offset(net));
    } else if (is_simple_snarl(record.get_record_type()) && is_chain(net)) {
        return SimpleSnarlRecord(net, &snarl_tree_records).get_node_is_reversed();
    } else {
        return record.get_is_reversed_in_parent();
    }
}
net_handle_t SnarlDistanceIndex::get_node_net_handle(const nid_t id, bool rev) const {
    RootRecord root_record (get_root(), &snarl_tree_records);
    size_t node_pointer_offset = get_node_pointer_offset(id, root_record.get_min_node_id(), root_record.get_connected_component_count());
    size_t record_offset = snarl_tree_records->at(node_pointer_offset);
    size_t node_record_offset = snarl_tree_records->at(node_pointer_offset+1);
    return get_net_handle_from_values(record_offset, rev ? END_START : START_END, NODE_HANDLE, node_record_offset); 
}

size_t SnarlDistanceIndex::get_max_tree_depth() const {
    return RootRecord(get_root(), &snarl_tree_records).get_max_tree_depth();
}

size_t SnarlDistanceIndex::get_depth(const net_handle_t& net) const {
    if (is_root(net)) {
        return 0;
    } else if (is_simple_snarl(SnarlTreeRecord(net, &snarl_tree_records).get_record_type())){
        //If this is a simple snarl, then it can be a node, snarl, or chain

        //The depth of the snarl's parent chain
        size_t snarl_depth = get_depth(get_parent(get_net_handle_from_values(get_record_offset(net), START_END, SNARL_HANDLE, 1)));
        if (is_snarl(net) || is_sentinel(net)) {
            //If this is a snarl or sentinel of the snarl, then it has the depth of its parent
            return snarl_depth;
        } else {
            //If it is the snarl's child node or chain, then add 1
            return snarl_depth + 1;
        }
    } else if (is_snarl(net) || is_sentinel(net) || is_node(net)) {
        return get_depth(get_parent(net));
    } else if (is_trivial_chain(net)) {
        return get_depth(get_parent(net)) + 1;
    } else if (is_chain(net)) {
        return ChainRecord(net, &snarl_tree_records).get_depth();
    } else {
        throw runtime_error("error: Unknown handle type");
    }
        
}

size_t SnarlDistanceIndex::get_connected_component_number(const net_handle_t& net) const {
    if (get_record_offset(net) == 0) {
        throw runtime_error("error: trying to get the connected component number of the root");
    }
    net_handle_t child = net;
    net_handle_t parent = is_root(net) ? net : get_parent(net);
    while (!is_root(parent)) {
        child = parent;
        parent = get_parent(child);
    }
    if (get_record_offset(parent) == 0) {
        //If the parent is actually the root
        return SnarlTreeRecord(child, &snarl_tree_records).get_rank_in_parent();
    } else {
        //Otherwise, it must be a root-level snarl pretending to be a root
        return SnarlTreeRecord(parent, &snarl_tree_records).get_rank_in_parent();
    }
}


net_handle_t SnarlDistanceIndex::get_handle_from_connected_component(size_t num) const {
    size_t child_offset = snarl_tree_records->at(ROOT_RECORD_SIZE + num);
    net_handle_record_t type = SnarlTreeRecord(child_offset, &snarl_tree_records).get_record_handle_type();
    if (type == NODE_HANDLE) {
        //If this child is a node, then pretend it's a chain
        return get_net_handle_from_values(child_offset, START_END, CHAIN_HANDLE);
    } else if (type == SNARL_HANDLE) {
        return get_net_handle_from_values(child_offset, START_END, ROOT_HANDLE);
    } else {
        return get_net_handle_from_values(child_offset, START_END, type);
    }
}



bool SnarlDistanceIndex::has_connectivity(const net_handle_t& net, endpoint_t start, endpoint_t end) const {
    SnarlTreeRecord record(net, &snarl_tree_records);
    return record.has_connectivity(start, end);
}

SnarlDistanceIndex::SnarlTreeRecord::SnarlTreeRecord (size_t pointer, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = pointer;
    records = tree_records;
 
}

SnarlDistanceIndex::SnarlTreeRecord::SnarlTreeRecord (const net_handle_t& net, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = get_record_offset(net);
    records = tree_records;
}


bool SnarlDistanceIndex::SnarlTreeRecord::has_connectivity(connectivity_t connectivity) const {
    if (connectivity == START_START) {
        return is_start_start_connected();
    } else if (connectivity == START_END || connectivity == END_START) {
        return is_start_end_connected();
    } else if (connectivity == START_TIP || connectivity == TIP_START) {
        return is_start_tip_connected();
    } else if (connectivity == END_END) {
        return is_end_end_connected();
    } else if (connectivity == END_TIP || connectivity == TIP_END) {
        return is_end_tip_connected();
    } else if (connectivity == TIP_TIP) {
        return is_tip_tip_connected();
    } else {
        throw runtime_error("error: Invalid connectivity");
    }
}

size_t SnarlDistanceIndex::get_prefix_sum_value(const net_handle_t& net) const {
    return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_prefix_sum(get_node_record_offset(net));
}
size_t SnarlDistanceIndex::get_max_prefix_sum_value(const net_handle_t& net) const {
    return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_max_prefix_sum(get_node_record_offset(net));
}

size_t SnarlDistanceIndex::get_forward_loop_value(const net_handle_t& net) const {
    return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_forward_loop(get_node_record_offset(net));
}

size_t SnarlDistanceIndex::get_reverse_loop_value(const net_handle_t& net) const {
    return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_reverse_loop(get_node_record_offset(net));
}

size_t SnarlDistanceIndex::get_chain_component(const net_handle_t& net, bool get_end) const {
    return TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_chain_component(get_node_record_offset(net), get_end);
}


size_t SnarlDistanceIndex::SnarlTreeRecord::get_min_length() const {
    record_t type = get_record_type();
    size_t val;
    if (is_any_root(type)) {
        throw runtime_error("error: trying to find the length of the root");
    } else if (!has_distances(type)) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    } else if (is_node(type)) {
        return (*records)->at(record_offset + NODE_LENGTH_OFFSET);
    } else if (is_trivial_snarl(type)) {
        size_t last_node_offset = TrivialSnarlRecord(record_offset, records).get_node_count()-1;
        return (*records)->at(record_offset + DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (last_node_offset*2) + 1);
    } else if (is_simple_snarl(type))  {
        size_t raw_val =  (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET);
        return (raw_val >> 11) & ((1 << 11) - 1);
    } else if (is_nonroot_nonsimple_snarl(type))  {
        val =  (*records)->at(record_offset + SNARL_MIN_LENGTH_OFFSET);
    } else if (is_chain(type))  {
        val = (*records)->at(record_offset + CHAIN_MIN_LENGTH_OFFSET);
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
    return val == 0 ? std::numeric_limits<size_t>::max() : val - 1;
};

size_t SnarlDistanceIndex::SnarlTreeRecord::get_max_length() const {
    record_t type = get_record_type();
    size_t val;
    if (is_any_root(type)) {
        throw runtime_error("error: trying to find the length of the root");
    } else if (!has_distances(type)) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    } else if (is_node(type)) {
        return (*records)->at(record_offset + NODE_LENGTH_OFFSET);
    } else if (is_trivial_snarl(type)) {
        size_t last_node_offset = TrivialSnarlRecord(record_offset, records).get_node_count()-1;
        return (*records)->at(record_offset + DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (last_node_offset*2) + 1);
    } else if (is_simple_snarl(type))  {
        size_t raw_val =  (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET);
        return raw_val & ((1 << 11) - 1);
    } else if (is_nonroot_nonsimple_snarl(type))  {
        val = (*records)->at(record_offset + SNARL_MAX_LENGTH_OFFSET);
    } else if (is_chain(type))  {
        val = (*records)->at(record_offset + CHAIN_MAX_LENGTH_OFFSET);
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }

    return val == 0 ? std::numeric_limits<size_t>::max() : val - 1;
};

size_t SnarlDistanceIndex::SnarlTreeRecord::get_rank_in_parent() const {
    record_t type = get_record_type();
    if (is_node(type)) {
        return (*records)->at(record_offset + NODE_RANK_OFFSET);
    } else if (is_trivial_snarl(type))  {
        throw runtime_error("error: node ranks need the node offsets");
    } else if (is_root_snarl(type)) {
        //For root snarls, the rank gets stored in the length slot
        return (*records)->at( record_offset + SNARL_MIN_LENGTH_OFFSET);
    } else if (is_nonroot_nontrivial_snarl(type))  {
        return record_offset;
    } else if (is_chain(type))  {
        return (*records)->at(record_offset + CHAIN_RANK_OFFSET) >> 1;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
};
bool SnarlDistanceIndex::SnarlTreeRecord::get_is_reversed_in_parent() const {
    record_t type = get_record_type();
    if (is_node(type)) {
        return false;
    } else if (is_trivial_snarl(type)) {
        cerr << "warning: Getting orientation of a trivial snarl" << endl;
        return false;
    } else if (is_nontrivial_snarl(type))  {
        return false;
    } else if (is_chain(type))  {
        return (*records)->at(record_offset + CHAIN_RANK_OFFSET) & 1;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
};

handlegraph::nid_t SnarlDistanceIndex::SnarlTreeRecord::get_start_id() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        throw runtime_error("error: trying to get the start node of the root");
    } else if (is_node(type)) {
        //cerr << "warning: Looking for the start of a node" << endl;
        return (*records)->at(record_offset + NODE_ID_OFFSET);
    } else if (is_trivial_snarl(type)) {
        TrivialSnarlRecord trivial_snarl_record(record_offset, records);
        return trivial_snarl_record.get_node_id(0);
    } else if (is_nonroot_nontrivial_snarl(type))  {
        //To get the start node of a snarl, get the thing to the left of it in the chain
        ChainRecord parent_record (get_parent_record_offset(), records);
        size_t node_offset = is_simple_snarl(type) ? 1 : 0;
        net_handle_t next_node_in_chain = parent_record.get_next_child(get_net_handle_from_values(record_offset, START_END, SNARL_HANDLE, node_offset), true);
        return TrivialSnarlRecord(get_record_offset(next_node_in_chain), records).get_node_id(get_node_record_offset(next_node_in_chain));
    } else if (is_chain(type))  {
        return ((*records)->at(record_offset + CHAIN_START_NODE_OFFSET)) >> 1;
    } else if (is_root_snarl(type)) {
        throw runtime_error("error: trying to find the start node of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
}
bool SnarlDistanceIndex::SnarlTreeRecord::get_start_orientation() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        throw runtime_error("error: trying to get the start node of the root");
    } else if (is_node(type) || is_trivial_snarl(type)) {
        //cerr << "warning: Looking for the start of a node" << endl;
        return false;
    } else if (is_nonroot_nontrivial_snarl(type))  {
        ChainRecord parent_record (get_parent_record_offset(), records);
        size_t node_offset = is_simple_snarl(type) ? 1 : 0;
        //Get the next node in the chain (going left)
        //The handle will be pointing in the direction we just moved, so if it is going START_END, then it is reversed
        net_handle_t next_node_in_chain = parent_record.get_next_child(get_net_handle_from_values(record_offset, START_END, SNARL_HANDLE, node_offset), true);
        return get_end_endpoint( get_connectivity(next_node_in_chain)) == END;
    } else if (is_chain(type))  {
        return ((*records)->at(record_offset + CHAIN_START_NODE_OFFSET)) & 1;
    } else if (is_root_snarl(type)) {
        throw runtime_error("error: trying to find the start node of a root snarl");
    }else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
}
handlegraph::nid_t SnarlDistanceIndex::SnarlTreeRecord::get_end_id() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        throw runtime_error("error: trying to get the end node of the root");
    } else if (is_node(type)) {
        //cerr << "warning: Looking for the end of a node" << endl;
        //Offset of the start of the node vector
        return (*records)->at(record_offset + NODE_ID_OFFSET);
    } else if (is_trivial_snarl(type)) {
        TrivialSnarlRecord trivial_snarl_record(record_offset, records);
        return trivial_snarl_record.get_node_id(trivial_snarl_record.get_node_count()-1);
    } else if (is_nonroot_nontrivial_snarl(type))  {
        //For a snarl, walk right in the chain
        ChainRecord parent_record (get_parent_record_offset(), records);
        size_t node_offset = is_simple_snarl(type) ? 1 : 0;
        net_handle_t next_node_in_chain = parent_record.get_next_child(get_net_handle_from_values(record_offset, START_END, SNARL_HANDLE, node_offset), false);
        return TrivialSnarlRecord(get_record_offset(next_node_in_chain), records).get_node_id(get_node_record_offset(next_node_in_chain));
    } else if (is_chain(type))  {
        return ((*records)->at(record_offset + CHAIN_END_NODE_OFFSET)) >> 1;
    } else if (is_root_snarl(type)) {
        throw runtime_error("error: trying to find the end node of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
}

bool SnarlDistanceIndex::SnarlTreeRecord::get_end_orientation() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        throw runtime_error("error: trying to get the end node of the root");
    } else if (is_node(type) || is_trivial_snarl(type)) {
        //cerr << "warning: Looking for the end of a node" << endl;
        //Offset of the start of the node vector
        return false;
    } else if (is_nonroot_nontrivial_snarl(type))  {
        ChainRecord parent_record (get_parent_record_offset(), records);
        size_t node_offset = is_simple_snarl(type) ? 1 : 0;
        //Get the next node in the chain (going right)
        net_handle_t next_node_in_chain = parent_record.get_next_child(get_net_handle_from_values(record_offset, START_END, SNARL_HANDLE, node_offset), false);
        //The handle will be pointing in the direction we just moved, so if it is going END_START, then it is reversed
        return get_end_endpoint( get_connectivity(next_node_in_chain)) == START;
    } else if (is_chain(type))  {
        return ((*records)->at(record_offset + CHAIN_END_NODE_OFFSET)) & 1;
    } else if (is_root_snarl(type)) {
        throw runtime_error("error: trying to find the end node of a root snarl");
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
}



bool SnarlDistanceIndex::SnarlTreeRecord::has_connectivity(endpoint_t start, endpoint_t end){
    return has_connectivity(endpoints_to_connectivity(start, end));
}

bool SnarlDistanceIndex::has_external_connectivity(const net_handle_t& net, endpoint_t start, endpoint_t end) const {
    size_t tag = snarl_tree_records->at(get_record_offset(net));
    if (start == START && end == START) {
        return is_externally_start_start_connected(tag);
    } else if (start == END && end == END) {
        return is_externally_end_end_connected(tag);
    } else if ( (start == START && end == END ) || (start == END && end == START)) {
        return is_externally_start_end_connected(tag);
    } else {
        return false;
    }
}

size_t SnarlDistanceIndex::SnarlTreeRecord::get_parent_record_offset() const {
    record_t type = get_record_type();
    if (type == ROOT) {
        return 0;
    } else if (is_node(type)) {
        return ((*records)->at(record_offset + NODE_PARENT_OFFSET));
    } else if (is_trivial_snarl(type)) {
        return (*records)->at(record_offset + TRIVIAL_SNARL_PARENT_OFFSET);
    } else if (is_simple_snarl(type)) {
        return (*records)->at(record_offset + SIMPLE_SNARL_PARENT_OFFSET);
    } else if (is_nonsimple_snarl(type))  {
        return ((*records)->at(record_offset + SNARL_PARENT_OFFSET));
    } else if (is_chain(type))  {
        return ((*records)->at(record_offset + CHAIN_PARENT_OFFSET));
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }
};

SnarlDistanceIndex::SnarlTreeRecordWriter::SnarlTreeRecordWriter (size_t pointer, bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = pointer;
    records = tree_records;
 
}

SnarlDistanceIndex::SnarlTreeRecordWriter::SnarlTreeRecordWriter (const net_handle_t& net, 
                                                          bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = get_record_offset(net);
    records = tree_records;
}




SnarlDistanceIndex::record_t SnarlDistanceIndex::SnarlTreeRecordWriter::get_record_type() const {
    return static_cast<record_t>((*records)->at(record_offset) >> 9);
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_start_start_connected() {

    (*records)->at(record_offset) = (*records)->at(record_offset) | 32;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_start_end_connected() {

    (*records)->at(record_offset) = (*records)->at(record_offset) | 16;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_start_tip_connected() {

    (*records)->at(record_offset) = (*records)->at(record_offset) | 8;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_end_end_connected() {

    (*records)->at(record_offset) = (*records)->at(record_offset) | 4;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_end_tip_connected() {

    (*records)->at(record_offset) = (*records)->at(record_offset) | 2;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_tip_tip_connected() {

    (*records)->at(record_offset) = (*records)->at(record_offset) | 1;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_externally_start_end_connected() {

    (*records)->at(record_offset) = (*records)->at(record_offset) | 64;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_externally_start_start_connected() const {
    (*records)->at(record_offset) = (*records)->at(record_offset) | 128;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_externally_end_end_connected() const {
    (*records)->at(record_offset) = (*records)->at(record_offset) | 256;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_record_type(record_t type) {
    assert((*records)->at(record_offset) == 0);
    (*records)->at(record_offset) = ((static_cast<size_t>(type) << 9) | ((*records)->at(record_offset) & 511));
}


void SnarlDistanceIndex::SnarlTreeRecordWriter::set_min_length(size_t length) {
    record_t type = get_record_type();
    size_t offset;
    if (!has_distances(type)) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    } else if (is_any_root(type)) {
        throw runtime_error("error: set the length of a root snarl");
    } else if (is_node(type)) {
        offset = record_offset + NODE_LENGTH_OFFSET;
    } else if (is_nonroot_nonsimple_snarl(type))  {
        offset = record_offset + SNARL_MIN_LENGTH_OFFSET;
    } else if (is_chain(type))  {
        offset = record_offset + CHAIN_MIN_LENGTH_OFFSET;
    } else if (is_simple_snarl(type)) {
        if (length > (1 << 11)-1) {
            throw runtime_error("error: node length is too large");
        }
        size_t old_val =  (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET);
        size_t new_val = old_val | (length << 11);
        (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET)= new_val;
         return;  
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }

    (*records)->at(offset) = (length == std::numeric_limits<size_t>::max() ? 0 : length + 1);
};
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_max_length(size_t length) {
    record_t type = get_record_type();
    size_t offset;
    if (!has_distances(type)) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    } else if (is_any_root(type)) {
        throw runtime_error("error: set the length of a root snarl");
    } else if (is_node(type)) {
        throw runtime_error("error: set the max length of a node");
    } else if (is_nonroot_nonsimple_snarl(type))  {
        offset = record_offset + SNARL_MAX_LENGTH_OFFSET;
    } else if (is_simple_snarl(type)) {
        if (length > (1 << 11)-1) {
            throw runtime_error("error: node length is too large");
        }
        size_t old_val =  (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET);
        size_t new_val = old_val | length;
        (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET) = new_val;
        return;
    } else if (is_chain(type))  {
        offset = record_offset + CHAIN_MAX_LENGTH_OFFSET;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }

    (*records)->at(offset) = (length == std::numeric_limits<size_t>::max() ? 0 : length + 1);
};

void SnarlDistanceIndex::SnarlTreeRecordWriter::set_rank_in_parent(size_t rank) {
    record_t type = get_record_type();
    size_t offset;
    if (is_node(type)) {
        offset = record_offset + NODE_RANK_OFFSET;
    } else if (is_root_snarl(type)) {
        (*records)->at(record_offset + SNARL_MIN_LENGTH_OFFSET) = rank;
        return;
    } else if (is_nonroot_nontrivial_snarl(type))  {
        cerr << "SETTING THE RANK OF A " << stringify(type) << " SNARL WHICH I'M PRETTY SURE DOESN'T MEAN ANYTHING" << endl;
        return;
    } else if (is_chain(type))  {
        offset = record_offset + CHAIN_RANK_OFFSET;
    } else {
        throw runtime_error("error: trying to set rank of a snarl tree node of the wrong type: " + stringify(type));
    }

    bool rev = (*records)->at(offset) & 1;
    (*records)->at(offset) = (rank << 1) | rev;
};
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_is_reversed_in_parent(bool rev) {
    record_t type = get_record_type();
    size_t offset;
    if (is_nontrivial_snarl(type))  {
        return;
    } else if (is_chain(type))  {
        offset = record_offset + CHAIN_RANK_OFFSET;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }

    (*records)->at(offset) =  (((*records)->at(offset)>>1)<<1) | rev;
};
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_parent_record_offset(size_t pointer){
    record_t type = get_record_type();
    size_t offset;
    if (is_node(type)) {
        offset = record_offset + NODE_PARENT_OFFSET;
    }  else if (is_trivial_snarl(type)) {
        offset = record_offset + TRIVIAL_SNARL_PARENT_OFFSET;
    } else if (is_nonsimple_snarl(type))  {

        offset = record_offset + SNARL_PARENT_OFFSET;
    } else if (is_simple_snarl(type)) {
        offset = record_offset + SIMPLE_SNARL_PARENT_OFFSET;
    } else if (is_chain(type))  {
        offset = record_offset + CHAIN_PARENT_OFFSET;
    } else {
        throw runtime_error("error: trying to access a snarl tree node of the wrong type");
    }

    (*records)->at(offset) = pointer;
};
//Rev is true if the node is traversed backwards to enter the snarl
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_start_node(handlegraph::nid_t id, bool rev) {
    record_t type = get_record_type();
    size_t offset;
    if (type == ROOT || is_node(type) || is_trivial_snarl(type)) {
        throw runtime_error("error: trying to set the start node id of a node or root");
    } else if (is_nonroot_nontrivial_snarl(type)) {
        throw runtime_error("error: trying to set the start node id of a snarl");
    } else if (is_chain(type)) {
        offset = record_offset + CHAIN_START_NODE_OFFSET;
    } else if (is_root_snarl(type)) {
        throw runtime_error("error: set the start node of a root snarl");
    } else {
        throw runtime_error("error: trying to set start on a snarl tree node of the wrong type: " + stringify(type));
    }

    (*records)->at(offset) = (id << 1) | rev;
}
void SnarlDistanceIndex::SnarlTreeRecordWriter::set_end_node(handlegraph::nid_t id, bool rev) const {
    record_t type = get_record_type();
    size_t offset;
    if (type == ROOT || is_node(type) || is_trivial_snarl(type)) {
        throw runtime_error("error: trying to set the node id of a node or root");
    } else if (is_nonroot_nontrivial_snarl(type)) {
        throw runtime_error("error: trying to set the end node id of a snarl");
    } else if (is_chain(type)) {
        offset = record_offset + CHAIN_END_NODE_OFFSET;
    } else if (is_root_snarl(type)) {
        throw runtime_error("error: set the end node of a root snarl");
    } else {
        throw runtime_error("error: trying to set end on a snarl tree node of the wrong type: " + stringify(type));
    }

    (*records)->at(offset) = (id << 1) | rev;
}

SnarlDistanceIndex::RootRecord::RootRecord (size_t pointer, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = pointer;
    records = tree_records;
}


SnarlDistanceIndex::RootRecord::RootRecord (net_handle_t net, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = get_record_offset(net);
    records = tree_records;
}



bool SnarlDistanceIndex::RootRecord::for_each_child(const std::function<bool(const handlegraph::net_handle_t&)>& iteratee) const {
    assert(record_offset == 0);
    size_t connected_component_count = get_connected_component_count();
    for (size_t i = 0 ; i < connected_component_count ; i++) {

        size_t child_offset = (*records)->at(record_offset + ROOT_RECORD_SIZE + i);
        net_handle_record_t type = SnarlTreeRecord(child_offset, records).get_record_handle_type();
        record_t record_type = SnarlTreeRecord(child_offset, records).get_record_type();


        if (is_root_snarl(record_type)) {
            //This is a bunch of root components that are connected, so go through each
            SnarlRecord snarl_record(child_offset, records);
            if (! snarl_record.for_each_child(iteratee)) {
                return false;
            }
        } else {
            //Otherwise, it is a separate connected component
            net_handle_t child_handle;
            if (type == NODE_HANDLE) {
                //If this child is a node, then pretend it's a chain
                child_handle = get_net_handle_from_values(child_offset, START_END, CHAIN_HANDLE);
            } else {
                child_handle =  get_net_handle_from_values(child_offset, START_END, type);
            }
            if (!iteratee(child_handle)) {
                return false;
            }
        }
    }
    return true;
}

SnarlDistanceIndex::RootRecordWriter::RootRecordWriter (size_t pointer, size_t connected_component_count, size_t node_count, 
        size_t max_tree_depth, handlegraph::nid_t min_node_id, bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* records){

    record_offset = pointer;
    records = records;
    SnarlTreeRecordWriter::record_offset = pointer;
    SnarlTreeRecordWriter::records = records;
    RootRecord::record_offset = pointer;
    RootRecord::records = records;
    //Allocate memory for the root vector and for all of the nodes
    (*records)->resize((*records)->size() + ROOT_RECORD_SIZE + connected_component_count + (node_count*2));
    set_record_type(ROOT);
    set_min_node_id(min_node_id);
    set_node_count(node_count);
    set_max_tree_depth(max_tree_depth);
    set_connected_component_count(connected_component_count);
    set_version_number();
}

void SnarlDistanceIndex::RootRecordWriter::set_version_number() {
    (*records)->at(record_offset+VERSION_NUMBER_OFFSET) = VERSION_NUMBER_SENTINEL ^ CURRENT_VERSION_NUMBER;
}

void SnarlDistanceIndex::RootRecordWriter::set_connected_component_count(size_t connected_component_count) {

    (*records)->at(record_offset+COMPONENT_COUNT_OFFSET)=connected_component_count;
}
void SnarlDistanceIndex::RootRecordWriter::set_node_count(size_t node_count) {

    (*records)->at(record_offset+NODE_COUNT_OFFSET)=node_count;
}
void SnarlDistanceIndex::RootRecordWriter::set_max_tree_depth(size_t tree_depth) {

    (*records)->at(record_offset+MAX_TREE_DEPTH_OFFSET)=tree_depth;
}
void SnarlDistanceIndex::RootRecordWriter::set_min_node_id(handlegraph::nid_t node_id) {

    (*records)->at(record_offset+MIN_NODE_ID_OFFSET)=node_id;
}
void SnarlDistanceIndex::RootRecordWriter::add_component(size_t index, size_t offset) {

    (*records)->at(record_offset+ROOT_RECORD_SIZE+index) = offset;

}

SnarlDistanceIndex::SnarlRecord::SnarlRecord (size_t pointer, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = pointer;
    records = tree_records;
}

SnarlDistanceIndex::SnarlRecord::SnarlRecord (net_handle_t net, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records) :
    SnarlRecord(get_record_offset(net), tree_records) {
    // Nothing to do!
}


size_t SnarlDistanceIndex::SnarlRecord::distance_vector_size(record_t type, size_t node_count) {
    if (!is_nonsimple_snarl(type)) {
        throw runtime_error("error: trying to get size of distance matrix for something other than a snarl that would have one");
    }
    if (has_distances(type)) {
        if (is_oversized_snarl(type)) {
            // For oversized snarls, hhl_size for the hub labeling data is accounted separately.
            return 0;
        } else {
            // Must be a distanced non-simple snarl or root snarl.
            
            //For a normal min distance snarl just the distances between internal node sides
            size_t node_side_count = node_count * 2;
            size_t vector_size =  (((node_side_count+1)*node_side_count) / 2);
            return vector_size;
        }
    } else {
        //For a non-distanced snarl, its just the record size and the pointers to children
        return 0;
    }
}

size_t SnarlDistanceIndex::SnarlRecord::record_size (record_t type, size_t node_count, size_t hhl_size) {
    if (is_oversized_snarl(type)) {
      // Oversized snarls need the fixed-size header and the packed hub label
      // vector.

      return SNARL_RECORD_SIZE + hhl_size;
    } else {
      // Normal snarl records need the fixed-size header and the distance matrix
      return SNARL_RECORD_SIZE + distance_vector_size(type, node_count);
    }
}

size_t SnarlDistanceIndex::SnarlRecord::get_distance_start_start() const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + SNARL_DISTANCE_START_START_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}
size_t SnarlDistanceIndex::SnarlRecord::get_distance_end_end() const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + SNARL_DISTANCE_END_END_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}

SnarlDistanceIndex::SnarlRecordWriter::SnarlRecordWriter (size_t node_count, bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* records, record_t type, size_t hhl_size){
    //Constructor for making a new record, including allocating memory.
    //Assumes that this is the latest record being made, so pointer will be the end of
    //the array and we need to allocate extra memory past it

    records = records;
    record_offset = (*records)->size();
    SnarlTreeRecordWriter::record_offset = (*records)->size();
    SnarlTreeRecordWriter::records = records;
    SnarlRecord::record_offset = (*records)->size();
    SnarlRecord::records = records;
    
    //hhl_size only used for oversized snarls
    size_t extra_size = record_size(type, node_count, hhl_size);
    (*records)->resize((*records)->size() + extra_size);
    set_node_count(node_count);
    set_record_type(type);
    
}

SnarlDistanceIndex::SnarlRecordWriter::SnarlRecordWriter(bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* records, size_t pointer) {
    //Make a constructor for a snarl record that has already been allocated.
    //For adding children to an existing snarl record
    record_offset = pointer;
    records = records;
    SnarlTreeRecordWriter::records = records;
    SnarlTreeRecordWriter::record_offset = pointer;
    SnarlRecord::records = records;
    SnarlRecord::record_offset = pointer;
}


void SnarlDistanceIndex::SnarlRecordWriter::set_distance_start_start(size_t value) {
    (*records)->at(record_offset + SNARL_DISTANCE_START_START_OFFSET) = (value == std::numeric_limits<size_t>::max() ? 0 : value + 1);
}
void SnarlDistanceIndex::SnarlRecordWriter::set_distance_end_end(size_t value) {
    (*records)->at(record_offset + SNARL_DISTANCE_END_END_OFFSET) = (value == std::numeric_limits<size_t>::max() ? 0 : value + 1);
}


size_t SnarlDistanceIndex::SnarlRecord::get_distance_vector_offset(size_t rank1, bool right_side1, size_t rank2,
        bool right_side2, size_t node_count, record_t type) {
    if (!has_distances(type)) {
        throw runtime_error("error: trying to access distance in a distanceless snarl tree");
    }
    if (is_oversized_snarl(type)) {
        throw runtime_error("error: trying to access distance matrix in an oversized snarl");
    }
    if (!is_nonsimple_snarl(type)) {
        throw runtime_error("error: trying to access distance matrix in something other than a snarl that would have one");
    }

    //how many node sides in this snarl
    size_t node_side_count = node_count * 2;

    if (!is_root_snarl(type)) {
        //For distances snarls, the ranks 0 and 1 are for the start and end nodes. The distance
        //matrix in the snarl record is only for distances between internal nodes, so 0 and 1 are
        //never stored and we decrement the ranks by 2
        rank1 -= 2;
        rank2 -= 2;
    }

    //make sure we're looking at the correct node side
    rank1 = rank1 * 2;
    if (right_side1) {
        rank1 += 1;
    }
    rank2 = rank2 * 2;
    if (right_side2) {
        rank2 += 1;
    }
    

    //reverse order of ranks if necessary
    if (rank1 > rank2) {
        size_t tmp = rank1;
        rank1 = rank2;
        rank2 = tmp;
    }

    //normal distance index
    size_t k = node_side_count-rank1;
    return (((node_side_count+1) * node_side_count)/2) - (((k+1)*k) / 2) + rank2 - rank1;
}

size_t SnarlDistanceIndex::SnarlRecord::get_distance_vector_offset(size_t rank1, bool right_side1,
        size_t rank2, bool right_side2) const {
    return get_distance_vector_offset(rank1, right_side1, rank2, right_side2,
        get_node_count(), get_record_type());
}


bool SnarlDistanceIndex::SnarlRecord::for_each_child(const std::function<bool(const net_handle_t&)>& iteratee) const {
    size_t child_count = get_node_count();
    size_t child_record_offset = get_child_record_pointer();
    for (size_t i = 0 ; i < child_count ; i++) {
        size_t child_offset =  (*records)->at(child_record_offset + i);
        net_handle_t child_handle =  get_net_handle_from_values (child_offset, START_END, CHAIN_HANDLE);
        bool result = iteratee(child_handle);
        if (result == false) {
            return false;
        }
    }
    return true;
}

void SnarlDistanceIndex::SnarlRecordWriter::set_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2, size_t distance) {

    //Don't save internal distances for oversized snarls
    if (is_oversized_snarl(get_record_type()) && !(rank1 == 0 || rank1 == 1 || rank2 == 0 || rank2 == 1)){
        return;
    }

    //Offset of this particular distance in the distance vector
    size_t distance_vector_offset = get_distance_vector_offset(rank1, right_side1, rank2, right_side2);
    assert(distance_vector_offset <= distance_vector_size(get_record_type(), get_node_count())); 

    //Value we actually want to save
    size_t val = distance == std::numeric_limits<size_t>::max() ? 0 : distance+1;

    (*records)->at(distance_vector_offset+record_offset+SNARL_RECORD_SIZE) = val;
}

size_t SnarlDistanceIndex::SnarlRecord::get_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    if (is_oversized_snarl(get_record_type())) {
        throw runtime_error("error: trying to distance from an oversized snarl");
    }

    //Offset of this particular distance in the distance vector
    size_t distance_vector_offset = get_distance_vector_offset(rank1, right_side1, rank2, right_side2);

    size_t val = (*records)->at(distance_vector_offset+record_offset+SNARL_RECORD_SIZE);

    return  val == 0 ? std::numeric_limits<size_t>::max() : val-1;

}
size_t SnarlDistanceIndex::SnarlRecord::get_node_count() const {
    return (*records)->at(record_offset + SNARL_NODE_COUNT_OFFSET);
}

//Node count is the number of nodes in the snarl, not including boundary nodes
void SnarlDistanceIndex::SnarlRecordWriter::set_node_count(size_t node_count) {

    (*records)->at(record_offset + SNARL_NODE_COUNT_OFFSET) = node_count;
}

void SnarlDistanceIndex::SnarlRecordWriter::set_hhl_entry(size_t index, size_t value) {
    (*records)->at(record_offset + SNARL_RECORD_SIZE + index) = value;
}


size_t SnarlDistanceIndex::SnarlRecord::get_child_record_pointer() const {
    return (*records)->at(record_offset+SNARL_CHILD_RECORD_OFFSET) ;
}

void SnarlDistanceIndex::SnarlRecordWriter::set_child_record_pointer(size_t pointer) {

    (*records)->at(record_offset+SNARL_CHILD_RECORD_OFFSET) = pointer;
}

//Add a reference to a child of this snarl. Assumes that the index is completed up
//to here
void SnarlDistanceIndex::SnarlRecordWriter::add_child(size_t pointer){
    size_t start_i = (*records)->size();
    (*records)->resize(start_i+1);
    (*records)->at(start_i) = pointer;


}

SnarlDistanceIndex::SimpleSnarlRecord::SimpleSnarlRecord (size_t pointer, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records, size_t node){
    record_offset = pointer;
    records = tree_records;
    node_rank = node;
}

SnarlDistanceIndex::SimpleSnarlRecord::SimpleSnarlRecord (net_handle_t net, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = get_record_offset(net);
    records = tree_records;
    node_rank = get_node_record_offset(net) == 1 ? std::numeric_limits<size_t>::max() : get_node_record_offset(net);
    
}



size_t SnarlDistanceIndex::SimpleSnarlRecord::record_size() {
    
    return SIMPLE_SNARL_RECORD_SIZE + get_node_count()*2;
}
size_t SnarlDistanceIndex::SimpleSnarlRecord::get_node_count() const {

    return (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET) >> 22;
}

SnarlDistanceIndex::SimpleSnarlRecordWriter::SimpleSnarlRecordWriter (size_t node_count, bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* records, record_t type){

    //Constructor for making a new record, including allocating memory.
    //Assumes that this is the latest record being made, so pointer will be the end of
    //the array and we need to allocate extra memory past it
    records = records;
    record_offset = (*records)->size();
    SnarlTreeRecordWriter::record_offset = (*records)->size();
    SnarlTreeRecordWriter::records = records;
    SimpleSnarlRecord::record_offset = (*records)->size();
    SimpleSnarlRecord::records = records;
    
    (*records)->resize((*records)->size() + SIMPLE_SNARL_RECORD_SIZE + 2*node_count);
    set_node_count(node_count);
    set_record_type(type);
    set_node_count(node_count);
    set_start_end_connected();
}

SnarlDistanceIndex::SimpleSnarlRecordWriter::SimpleSnarlRecordWriter(bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* records, size_t pointer) {
    //Make a constructor for a snarl record that has already been allocated.
    //For adding children to an existing snarl record
    record_offset = pointer;
    records = records;
    SnarlTreeRecordWriter::record_offset = pointer;
    SnarlTreeRecordWriter::records = records;
    SimpleSnarlRecord::record_offset = pointer;
    SimpleSnarlRecord::records = records;

}


void SnarlDistanceIndex::SimpleSnarlRecordWriter::set_node_count(size_t node_count) {
    size_t old_value = (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET);
    size_t new_value = old_value | (node_count << 22);
    (*records)->at(record_offset + SIMPLE_SNARL_NODE_COUNT_AND_LENGTHS_OFFSET) = new_value;
}
nid_t SnarlDistanceIndex::SimpleSnarlRecord::get_node_id(size_t rank) const {
    rank = rank == std::numeric_limits<size_t>::max() ? node_rank : rank;
    assert(rank >=2);
    assert(rank != std::numeric_limits<size_t>::max());
    return (*records)->at(record_offset + SIMPLE_SNARL_RECORD_SIZE + ((rank-2)*2));
}
size_t SnarlDistanceIndex::SimpleSnarlRecord::get_node_length(size_t rank) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    rank = rank == std::numeric_limits<size_t>::max() ? node_rank : rank;
    assert(rank >=2);
    assert(rank != std::numeric_limits<size_t>::max());
    return (*records)->at(record_offset + SIMPLE_SNARL_RECORD_SIZE + ((rank-2)*2) + 1) >> 1;
}
bool SnarlDistanceIndex::SimpleSnarlRecord::get_node_is_reversed(size_t rank) const {
    rank = rank == std::numeric_limits<size_t>::max() ? node_rank : rank;
    assert(rank>=2);
    assert(rank != std::numeric_limits<size_t>::max());
    return (*records)->at(record_offset + SIMPLE_SNARL_RECORD_SIZE + ((rank-2)*2) + 1) & 1;
}
size_t SnarlDistanceIndex::SimpleSnarlRecord::get_distance(size_t rank1, bool right_side1, size_t rank2, bool right_side2) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    if (rank1 > rank2) {
        //Order the nodes
        size_t tmp = rank1; bool tmp_rev = right_side1;
        rank1 = rank2; right_side1 = right_side2;
        rank2 = tmp; right_side2 = tmp_rev;
    }
    if (rank1 == 0 && rank2 == 1) {
        return get_min_length();
    } else if (rank1 == rank2) {
        //There are no loops
        return std::numeric_limits<size_t>::max();
    } else if (rank1 == 0) {
        //If the first node is the start node, then the distance is 0 if the node is not reversed 
        //and we are looking at the left side (or is reversed and right side), inf otherwise 
        if( get_node_is_reversed(rank2) == right_side2) {
            return 0;
        } else {
            return std::numeric_limits<size_t>::max();
        }
    } else if (rank1 == 1) {
        //If the first node is the end node, then the opposite of above
        if( get_node_is_reversed(rank2) == right_side2) {
            return std::numeric_limits<size_t>::max();
        } else {
            return 0;
        }
    } else {
        //Otherwise there is no path between them
        return std::numeric_limits<size_t>::max();
    }
}
net_handle_t SnarlDistanceIndex::SimpleSnarlRecord::get_child_from_rank(const size_t& rank) const {
    return get_net_handle_from_values(record_offset, START_END, CHAIN_HANDLE, rank);
}
bool SnarlDistanceIndex::SimpleSnarlRecord::for_each_child(const std::function<bool(const net_handle_t&)>& iteratee) const {
    size_t node_count = get_node_count();
    for (size_t i = 0 ; i < node_count ; i++) {
        net_handle_t child_handle = get_net_handle_from_values(record_offset, START_END, CHAIN_HANDLE, i+2);
        bool result = iteratee(child_handle);
        if (result == false) {
            return false;
        }
    }
    return true;
}
void SnarlDistanceIndex::SimpleSnarlRecordWriter::add_child(size_t rank, nid_t node_id, size_t node_length, bool is_reversed) {
    //Let the root know where to find the node
    (*records)->at(get_node_pointer_offset(node_id,
                                          (*records)->at(MIN_NODE_ID_OFFSET),
                                          (*records)->at(COMPONENT_COUNT_OFFSET))) 
                = record_offset;
    (*records)->at(get_node_pointer_offset(node_id,
                                          (*records)->at(MIN_NODE_ID_OFFSET),
                                          (*records)->at(COMPONENT_COUNT_OFFSET))+1) 
                = rank;
    //The rank gets decremented because the first 2 ranks are always the start and end nodes
    rank = rank-2;
    (*records)->at(record_offset+SIMPLE_SNARL_RECORD_SIZE+(rank*2)) = node_id;
    (*records)->at(record_offset+SIMPLE_SNARL_RECORD_SIZE+(rank*2)+1) = (node_length << 1) | is_reversed;

}

SnarlDistanceIndex::NodeRecord::NodeRecord (size_t pointer, size_t node_offset, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    record_offset = pointer;
    records = tree_records;

}

SnarlDistanceIndex::NodeRecord::NodeRecord (net_handle_t net, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    records = tree_records;
    record_offset = get_record_offset(net);
    
}


nid_t SnarlDistanceIndex::NodeRecord::get_node_id() const {
    return (*records)->at(record_offset + NODE_ID_OFFSET);
}
size_t SnarlDistanceIndex::NodeRecord::get_node_length() const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    return get_min_length();
    //if (get_record_type()== TRIVIAL_SNARL) {
    //    if (node_record_offset == 1) {
    //        return (*records)->at(record_offset + TRIVIAL_SNARL_RECORD_SIZE+1);
    //    } else {
    //        return ((*records)->at(record_offset + TRIVIAL_SNARL_RECORD_SIZE + ((node_record_offset-1)*2) + 1) - 
    //                (*records)->at(record_offset + TRIVIAL_SNARL_RECORD_SIZE + ((node_record_offset-1)*2) - 1));
    //    }
    //} else {
    //    return (*records)->at(record_offset + NODE_LENGTH_OFFSET);
    //}
}
size_t SnarlDistanceIndex::NodeRecord::get_distance_left_start() {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + NODE_DISTANCE_LEFT_START_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}
size_t SnarlDistanceIndex::NodeRecord::get_distance_right_start() {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + NODE_DISTANCE_RIGHT_START_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}
size_t SnarlDistanceIndex::NodeRecord::get_distance_left_end() {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + NODE_DISTANCE_LEFT_END_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}
size_t SnarlDistanceIndex::NodeRecord::get_distance_right_end() {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + NODE_DISTANCE_RIGHT_END_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}

SnarlDistanceIndex::TrivialSnarlRecord::TrivialSnarlRecord (size_t offset, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){
    records = tree_records;
    record_offset = offset;
    
}


size_t SnarlDistanceIndex::TrivialSnarlRecord::get_node_count() const {
    return (*records)->at(record_offset+TRIVIAL_SNARL_NODE_COUNT_OFFSET);
}

tuple<size_t, size_t, size_t, size_t> SnarlDistanceIndex::TrivialSnarlRecord::get_chain_values(size_t node_rank) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t prefix_sum = (*records)->at(record_offset+TRIVIAL_SNARL_PREFIX_SUM_OFFSET);  
    size_t forward_loop = (*records)->at(record_offset+TRIVIAL_SNARL_FORWARD_LOOP_OFFSET);  
    size_t reverse_loop = (*records)->at(record_offset+TRIVIAL_SNARL_REVERSE_LOOP_OFFSET);  
    size_t component = (node_rank == 0 && record_offset == ChainRecord(get_parent_record_offset(), records).get_first_node_offset()) ? 0: (*records)->at(record_offset+TRIVIAL_SNARL_COMPONENT_OFFSET); 

    //Get the right prefix sum value
    prefix_sum = prefix_sum == 0 ? std::numeric_limits<size_t>::max() : prefix_sum - 1;
    if (node_rank != 0) {
        prefix_sum =  sum(prefix_sum,  
               (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE+ (node_rank*2)-1));
    } 

    /*Get the forward loop value
     */

    forward_loop = forward_loop == 0 ? std::numeric_limits<size_t>::max() : forward_loop - 1;
        
    //How many nodes in the snarl
    size_t last_node = get_node_count()-1; 
    //The right prefix sum of the last node, which is the length of the snarl
    size_t snarl_length = (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (last_node*2) + 1); 

    //The distance from the right side of this node to the end of the trivial snarl
    //the length of the snarl - the right prefix sum of the node
    size_t right_offset = snarl_length - (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (node_rank*2) + 1);

    forward_loop = sum(forward_loop, right_offset*2);

    /* Get the reverse loop value
     */
    reverse_loop = reverse_loop == 0 ? std::numeric_limits<size_t>::max() : reverse_loop - 1;
        
    //The offset of the left side of this node in the trivial snarl
    size_t left_offset = node_rank == 0 ? 0 : 
                          (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (node_rank*2) - 1);

    reverse_loop = sum(reverse_loop, left_offset*2);



    return make_tuple(prefix_sum, forward_loop, reverse_loop, component);

}
size_t SnarlDistanceIndex::TrivialSnarlRecord::get_max_prefix_sum(size_t node_rank) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t prefix_sum = (*records)->at(record_offset+TRIVIAL_SNARL_MAX_PREFIX_SUM_OFFSET);

    //Get the right prefix sum value
    prefix_sum = prefix_sum == 0 ? std::numeric_limits<size_t>::max() : prefix_sum - 1;
    if (node_rank != 0) {
        prefix_sum =  sum(prefix_sum,  
               (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE+ (node_rank*2)-1));
    } 

    return prefix_sum;

}
size_t SnarlDistanceIndex::TrivialSnarlRecord::get_prefix_sum(size_t node_rank) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t prefix_sum = (*records)->at(record_offset+TRIVIAL_SNARL_PREFIX_SUM_OFFSET);  
    prefix_sum = prefix_sum == 0 ? std::numeric_limits<size_t>::max() : prefix_sum - 1;
    if (node_rank != 0) {
        return sum(prefix_sum,  
               (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE+ (node_rank*2)-1));
    } else {
        return prefix_sum;
    }
}
size_t SnarlDistanceIndex::TrivialSnarlRecord::get_forward_loop(size_t node_rank) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t forward_loop = (*records)->at(record_offset+TRIVIAL_SNARL_FORWARD_LOOP_OFFSET);  
    forward_loop = forward_loop == 0 ? std::numeric_limits<size_t>::max() : forward_loop - 1;
        
    //How many nodes in the snarl
    size_t last_node = get_node_count()-1; 
    //The right prefix sum of the last node, which is the length of the snarl
    size_t snarl_length = (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (last_node*2) + 1); 

    //The distance from the right side of this node to the end of the trivial snarl
    //the length of the snarl - the right prefix sum of the node
    size_t right_offset = snarl_length - (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (node_rank*2) + 1);

    return sum(forward_loop, right_offset*2);
}
size_t SnarlDistanceIndex::TrivialSnarlRecord::get_reverse_loop(size_t node_rank) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t reverse_loop = (*records)->at(record_offset+TRIVIAL_SNARL_REVERSE_LOOP_OFFSET);  
    reverse_loop = reverse_loop == 0 ? std::numeric_limits<size_t>::max() : reverse_loop - 1;
        
    //The offset of the left side of this node in the trivial snarl
    size_t left_offset = node_rank == 0 ? 0 : 
                          (*records)->at(record_offset+DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (node_rank*2) - 1);

    return sum(reverse_loop, left_offset*2);
}
size_t SnarlDistanceIndex::TrivialSnarlRecord::get_chain_component(size_t node_rank, bool get_end) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    if (!get_end && node_rank == 0 && record_offset == ChainRecord(get_parent_record_offset(), records).get_first_node_offset()){
        //If we want the component of the first node in the chain
        return 0;
    }
    return (*records)->at(record_offset + TRIVIAL_SNARL_COMPONENT_OFFSET);
}

size_t SnarlDistanceIndex::TrivialSnarlRecord::get_node_length(size_t node_rank) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    if (node_rank == 0) {
        return (*records)->at(record_offset + DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (2*node_rank)+1);
    } else {
        return (*records)->at(record_offset + DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (2*node_rank)+1) - 
               (*records)->at(record_offset + DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (2*node_rank)-1);
    }
}
size_t SnarlDistanceIndex::TrivialSnarlRecord::get_rank_in_parent(size_t node_rank) const {
    return record_offset+node_rank;
}
nid_t SnarlDistanceIndex::TrivialSnarlRecord::get_node_id(size_t node_rank) const {
    return get_record_type() == TRIVIAL_SNARL
          ? (*records)->at(record_offset + DISTANCELESS_TRIVIAL_SNARL_RECORD_SIZE + node_rank) >> 1
          : (*records)->at(record_offset + DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (2*node_rank)) >> 1;
}
bool SnarlDistanceIndex::TrivialSnarlRecord::get_is_reversed_in_parent(size_t node_rank) const {
    return get_record_type() == TRIVIAL_SNARL
          ? (*records)->at(record_offset + DISTANCELESS_TRIVIAL_SNARL_RECORD_SIZE + node_rank) & 1
          : (*records)->at(record_offset + DISTANCED_TRIVIAL_SNARL_RECORD_SIZE + (2*node_rank)) & 1;
}

SnarlDistanceIndex::NodeRecordWriter::NodeRecordWriter (size_t pointer, size_t node_offset, record_t type,
             bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* snarl_records, nid_t node_id){
    records = snarl_records;
    record_offset =  pointer;
    SnarlTreeRecordWriter::record_offset = pointer;
    SnarlTreeRecordWriter::records = records;
    NodeRecord::record_offset = pointer;
    NodeRecord::records = records;
    
    (*records)->resize((*records)->size() + NODE_RECORD_SIZE); 
    set_record_type(type);
    set_start_end_connected();
    
    //Set the pointer for the node to this record
    //Tell the root where to find this node, and a 0 for the offset in a trivial snarl
    (*records)->at(get_node_pointer_offset(node_id,
                                           (*records)->at(MIN_NODE_ID_OFFSET),
                                           (*records)->at(COMPONENT_COUNT_OFFSET))) = pointer;
    (*records)->at(get_node_pointer_offset(node_id,
                                           (*records)->at(MIN_NODE_ID_OFFSET),
                                           (*records)->at(COMPONENT_COUNT_OFFSET))+1) = 0;

}


void SnarlDistanceIndex::NodeRecordWriter::set_node_id(nid_t value) {

    (*records)->at(record_offset + NODE_ID_OFFSET) = value;
}
void SnarlDistanceIndex::NodeRecordWriter::set_node_length(size_t length) {
    (*records)->at(record_offset + NODE_LENGTH_OFFSET) = length;
}
void SnarlDistanceIndex::NodeRecordWriter::set_rank_in_parent(size_t value) {
    (*records)->at(record_offset + NODE_RANK_OFFSET) = value;
}
void SnarlDistanceIndex::NodeRecordWriter::set_distance_left_start(size_t distance) {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + NODE_DISTANCE_LEFT_START_OFFSET) = (distance == std::numeric_limits<size_t>::max() ? 0 : distance + 1);
}
void SnarlDistanceIndex::NodeRecordWriter::set_distance_right_start(size_t distance) {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + NODE_DISTANCE_RIGHT_START_OFFSET) = (distance == std::numeric_limits<size_t>::max() ? 0 : distance + 1);
}
void SnarlDistanceIndex::NodeRecordWriter::set_distance_left_end(size_t distance) {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + NODE_DISTANCE_LEFT_END_OFFSET) = (distance == std::numeric_limits<size_t>::max() ? 0 : distance + 1);
}
void SnarlDistanceIndex::NodeRecordWriter::set_distance_right_end(size_t distance) {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + NODE_DISTANCE_RIGHT_END_OFFSET) = (distance == std::numeric_limits<size_t>::max() ? 0 : distance + 1);
}

SnarlDistanceIndex::TrivialSnarlRecordWriter::TrivialSnarlRecordWriter (size_t pointer, record_t type,
             bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* snarl_records, bool new_record){

    records = snarl_records;
    record_offset =  pointer;
    SnarlTreeRecordWriter::record_offset = pointer;
    SnarlTreeRecordWriter::records = records;
    TrivialSnarlRecord::record_offset = pointer;
    TrivialSnarlRecord::records = records;
    
    assert (is_trivial_snarl(type));
    
 
    if (new_record) {
        if (type == TRIVIAL_SNARL) {
            (*records)->resize(record_offset + DISTANCELESS_TRIVIAL_SNARL_RECORD_SIZE);
        } else {
            (*records)->resize(record_offset + DISTANCED_TRIVIAL_SNARL_RECORD_SIZE);
        }
        set_record_type(type);
        set_start_end_connected();
    }
}


void SnarlDistanceIndex::TrivialSnarlRecordWriter::set_node_count(size_t value) const {
    (*records)->at(record_offset+TRIVIAL_SNARL_NODE_COUNT_OFFSET)=value;
}
void SnarlDistanceIndex::TrivialSnarlRecordWriter::set_prefix_sum(size_t value) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    value = value == std::numeric_limits<size_t>::max() ? 0 : value +1;
    (*records)->at(record_offset + TRIVIAL_SNARL_PREFIX_SUM_OFFSET) = value;
}
void SnarlDistanceIndex::TrivialSnarlRecordWriter::set_max_prefix_sum(size_t value) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    value = value == std::numeric_limits<size_t>::max() ? 0 : value +1;
    (*records)->at(record_offset + TRIVIAL_SNARL_MAX_PREFIX_SUM_OFFSET) = value;
}
void SnarlDistanceIndex::TrivialSnarlRecordWriter::set_forward_loop(size_t value) const { 
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    value = value == std::numeric_limits<size_t>::max() ? 0 : value +1;
    (*records)->at(record_offset + TRIVIAL_SNARL_FORWARD_LOOP_OFFSET) = value;
}
void SnarlDistanceIndex::TrivialSnarlRecordWriter::set_reverse_loop(size_t value) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    value = value == std::numeric_limits<size_t>::max() ? 0 : value +1;
    //Only do this for the first node
    (*records)->at(record_offset + TRIVIAL_SNARL_REVERSE_LOOP_OFFSET) = value;
    
}
void SnarlDistanceIndex::TrivialSnarlRecordWriter::set_chain_component(size_t value) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + TRIVIAL_SNARL_COMPONENT_OFFSET) = value;
}

SnarlDistanceIndex::ChainRecord::ChainRecord (size_t pointer, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records){

    record_offset = pointer;
    records = tree_records;

}       

SnarlDistanceIndex::ChainRecord::ChainRecord (net_handle_t net, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records) :
    ChainRecord(get_record_offset(net), tree_records) {
    // Nothing to do!
}   

SnarlDistanceIndex::ChainRecord::ChainRecord (net_handle_t net, const bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* tree_records, size_t tag){
    record_offset = get_record_offset(net);
    records = tree_records;
    
}   



size_t SnarlDistanceIndex::ChainRecord::get_node_count() const {
    if (get_record_handle_type() == NODE_HANDLE) {
        return 1;
    } else {
        return (*records)->at(record_offset + CHAIN_NODE_COUNT_OFFSET);
    }
}

tuple<size_t, bool, size_t> SnarlDistanceIndex::ChainRecord::get_last_child_offset() const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get children of a node");
    } else {
        size_t val = (*records)->at(record_offset + CHAIN_LAST_CHILD_OFFSET);
        size_t record_offset = val>>2;
        bool is_snarl = (val>>1) & 1;
        size_t node_record_offset = is_snarl ? 0 : TrivialSnarlRecord(record_offset, records).get_node_count()-1;
        return make_tuple(record_offset, is_snarl, node_record_offset);
    }
}
bool SnarlDistanceIndex::ChainRecord::get_is_looping_chain_connected_backwards() const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get children of a node");
    } else {
        size_t val = (*records)->at(record_offset + CHAIN_LAST_CHILD_OFFSET);
        return val & 1;
    }
}


size_t SnarlDistanceIndex::ChainRecord::get_distance(size_t rank1, bool left_side1, size_t node_length1, 
    size_t prefix_sum1, size_t forward_loop1, size_t reverse_loop1, size_t component1, size_t end_component1,
    size_t rank2, bool left_side2, size_t node_length2,
    size_t prefix_sum2, size_t forward_loop2, size_t reverse_loop2, size_t component2, size_t end_component2) const { 
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }

    //If 1 comes after 2, swap them
    if (rank1 > rank2) {
        size_t tmp_rank=rank1; bool tmp_side=left_side1; size_t tmp_len=node_length1;
        size_t tmp_pre=prefix_sum1; size_t tmp_fd = forward_loop1; size_t tmp_rev = reverse_loop1; size_t tmp_comp = component1;
        size_t tmp_end_comp = end_component1;
        rank1 = rank2; left_side1 = left_side2; node_length1 = node_length2;
        prefix_sum1 = prefix_sum2; forward_loop1 = forward_loop2; reverse_loop1 = reverse_loop2; component1 = component2; end_component1 = end_component2;
        rank2 = tmp_rank; left_side2 = tmp_side; node_length2 = tmp_len;
        prefix_sum2 = tmp_pre; forward_loop2 = tmp_fd; reverse_loop2 = tmp_rev; component2 = tmp_comp; end_component2 = tmp_end_comp;
    } 
    bool is_looping_chain = get_start_id() == get_end_id();
    
    record_t record_type = get_record_type();
    if (record_type == MULTICOMPONENT_CHAIN) {
        if (component1 != component2) {
            if (is_looping_chain) {

                //If this is a looping chain, then the first/last node could be in two
                //components
                return get_distance_taking_chain_loop(rank1, left_side1, node_length1, 
                            prefix_sum1, forward_loop1, reverse_loop1, end_component1,
                            rank2, left_side2, node_length2,
                            prefix_sum2, forward_loop2, reverse_loop2, end_component2);
            } else {
                return std::numeric_limits<size_t>::max();
            }
        }
    } else if (!is_chain(record_type)) {
        std::cerr << "Warning: weird record type for chain: " << stringify(record_type) << std::endl;
    }


    size_t distance;

    if (!left_side1 && left_side2) {
        //Right of 1 and left of 2, so a simple forward traversal of the chain
        if (rank1 == rank2) {
            //If these are the same node, then the path would need to go around the node
            distance = sum(sum(forward_loop1,reverse_loop2),node_length1);
        } else {
            distance = minus(prefix_sum2 - prefix_sum1, node_length1);
        }
    } else if (!left_side1 && !left_side2) {
        //Right side of 1 and right side of 2
        if (rank1 == rank2) {
            distance = forward_loop2;
        } else {
            distance = minus( sum(sum(prefix_sum2 - prefix_sum1, node_length2), forward_loop2),
                         node_length1);
        }
    } else if (left_side1 && left_side2) {
        //Left side of 1 and left side of 2
        if (rank1 == rank2) {
            distance = reverse_loop1;

        } else {
            distance = sum(prefix_sum2 - prefix_sum1, reverse_loop1);
        }
    } else {
        //Left side of 1 and right side of 2
        distance = sum(sum(sum(prefix_sum2 - prefix_sum1, reverse_loop1),
                        forward_loop2), node_length2);

    }
    if (is_looping_chain) {
        distance = std::min(distance, get_distance_taking_chain_loop(rank1, left_side1, node_length1, 
                            prefix_sum1, forward_loop1, reverse_loop1, end_component1,
                            rank2, left_side2, node_length2,
                            prefix_sum2, forward_loop2, reverse_loop2, end_component2));
    }
    return distance;
}


size_t SnarlDistanceIndex::ChainRecord::get_distance_taking_chain_loop(size_t rank1, bool left_side1, size_t node_length1, 
    size_t prefix_sum1, size_t forward_loop1, size_t reverse_loop1, size_t component1,
    size_t rank2, bool left_side2, size_t node_length2,
    size_t prefix_sum2, size_t forward_loop2, size_t reverse_loop2, size_t component2) const {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    //This is only called by get_distance, so the nodes should be ordered

    /*Note: Because we assume that the nodes are ordered and that we want to take the loop in the chain,
     * we leave the start node going left (either by going left or taking a loop to the right), and
     * enter the end node from the right (from the right or from the left and then looping)
     */
     if (get_record_type() == MULTICOMPONENT_CHAIN) {
        size_t last_component = TrivialSnarlRecord(get_first_node_offset() , records).get_chain_component(0, true);
        bool first_in_first_component = component1 == 0 || component1 == last_component;
        bool second_in_first_component = component2 == 0 || component2 == last_component;
        bool can_loop = get_is_looping_chain_connected_backwards();

        if (!can_loop || !first_in_first_component || ! second_in_first_component) {
            //If this is a multicomponent chain then it can only reach around backwards if both nodes
            //are in the first (last) component
            return std::numeric_limits<size_t>::max();
        } 

        //The first component of a multicomponent chain is the same as the last component
        //If we want to find the distance ending at the right side of a node in the first component
        //or starting at the left side of a node in the last component, then the chain has to be connected.
        //This is because this function assumes that we're taking the loop around the back of the chain,
        //and taking one of these paths would require the chain to be connected. 
        //eg. if we want to end at the right side of a node in the first component, then the start node must have
        //also been in the first component, so we'd go left, take the chain loop, then walk through the chain 
        //backwards from the end boundary node of the chain to the end node. If the end node is in the first component,
        //then the chain must be connected for the path to be valid
        bool end_at_right_of_first = component2 == 0 && !left_side2;
        bool start_at_left_of_last = component1 == last_component && left_side1; 
        if ((end_at_right_of_first || start_at_left_of_last)) {
            return std::numeric_limits<size_t>::max();
        }
    }


    size_t distance;

    if (!left_side1 && left_side2) {
        //Right of 1 and left of 2, so a simple forward traversal of the chain
        //loop forward from the first node, from the start of the chain to the first
        //node, from the end of the node to the second node, and the reverse loop of the second
        distance = sum(sum(sum(sum(forward_loop1, node_length1), prefix_sum1),
                        minus(get_min_length(), prefix_sum2)),
                        reverse_loop2);
    } else if (!left_side1 && !left_side2) {
        //Right side of 1 and right side of 2

        //Check distance for taking loop in chain: loop forward from the first node, from the start of the
        //chain to the first node, from the end of the node to the second node
        distance = sum(sum(sum(forward_loop1, node_length1), prefix_sum1),
                        minus(minus(get_min_length(), prefix_sum2), node_length2));
    } else if (left_side1 && left_side2) {
        //Left side of 1 and left side of 2

        //from the first node left to the start, around the
        //chain loop, then the reverse loop of the second node
        //This assumes that the length of the chain only includes the start/end node's length once,
        //which it does but might change
        distance = sum(sum(prefix_sum1,
                        minus(get_min_length(), prefix_sum2)),reverse_loop2);
    } else {
        //Left side of 1 and right side of 2

        //Check the distance going backwards around the chain
        distance = sum(prefix_sum1,
                        minus(minus(get_min_length(), prefix_sum2),node_length2));
    }
    return distance;
}

size_t SnarlDistanceIndex::ChainRecord::get_distance_left_start() {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + CHAIN_DISTANCE_LEFT_START_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}
size_t SnarlDistanceIndex::ChainRecord::get_distance_right_start() {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + CHAIN_DISTANCE_RIGHT_START_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}
size_t SnarlDistanceIndex::ChainRecord::get_distance_left_end() {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + CHAIN_DISTANCE_LEFT_END_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}
size_t SnarlDistanceIndex::ChainRecord::get_distance_right_end() {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    size_t stored_value = (*records)->at(record_offset + CHAIN_DISTANCE_RIGHT_END_OFFSET);
    return stored_value == 0 ? std::numeric_limits<size_t>::max() : stored_value - 1;
}

size_t SnarlDistanceIndex::ChainRecord::get_first_node_offset() const {
    if (get_record_handle_type() == NODE_HANDLE) {
        throw runtime_error("error: Trying to get chain traversal from a node");
    }
    return record_offset + CHAIN_RECORD_SIZE + 1;
}

net_handle_t SnarlDistanceIndex::ChainRecord::get_next_child(const net_handle_t& net_handle, bool go_left) const {
    //get the next child in the chain. net_handle must point to a snarl or node in the chain

    if (get_record_handle_type() == NODE_HANDLE) {
        //If this "chain" record is a node pretending to be a chain, then there is no next child
        return net_handle;
    }


    //If we are inside a trivial snarl and we aren't at the ends, then we just walk through the trivial snarl
    if (get_handle_type(net_handle) == NODE_HANDLE) {
        //If this is a node in a trivial snarl
        if (go_left && get_node_record_offset(net_handle) != 0) {
            //If we are going left and this is not the first node in the trivial snarl,
            //then keep everything the same but decrement the node record offset by one
            //to move one node left in the trivial snarl
            connectivity_t connectivity = get_connectivity(net_handle);
            if (get_node_record_offset(net_handle) == 1 && get_record_offset(net_handle) == get_first_node_offset() &&
                get_start_id() == get_end_id()) {
                //If this is the second thing in a looping chain and we're walking left to the first thing in the chain
                //then remember that the next handle is walking out of the chain
                connectivity = endpoints_to_connectivity(get_end_endpoint(connectivity), 
                                                                get_start_endpoint(connectivity));
            }
            if (TrivialSnarlRecord(get_record_offset(net_handle), records).get_is_reversed_in_parent(get_node_record_offset(net_handle)) !=
                TrivialSnarlRecord(get_record_offset(net_handle), records).get_is_reversed_in_parent(get_node_record_offset(net_handle)-1)) {
                //If this node is facing the opposite direction of the next nod in the chain
                connectivity = endpoints_to_connectivity(get_end_endpoint(connectivity), 
                                                                get_start_endpoint(connectivity));

            }
            return get_net_handle_from_values(get_record_offset(net_handle), connectivity, 
                                  get_handle_type(net_handle), get_node_record_offset(net_handle)-1);

        } else if (!go_left && get_node_record_offset(net_handle) != 
                      TrivialSnarlRecord(get_record_offset(net_handle), records).get_node_count()-1){

            //If we are going right and this is not the last node in the trivial snarl
            //then keep everything the same but increment the node record offset by one
            //to move one node right in the trivial snarl
            connectivity_t connectivity = get_connectivity(net_handle);
            if (TrivialSnarlRecord(get_record_offset(net_handle), records).get_is_reversed_in_parent(get_node_record_offset(net_handle)) !=
                TrivialSnarlRecord(get_record_offset(net_handle), records).get_is_reversed_in_parent(get_node_record_offset(net_handle)+1)) {
                //If this node is facing the opposite direction of the next nod in the chain
                connectivity = endpoints_to_connectivity(get_end_endpoint(connectivity), 
                                                                get_start_endpoint(connectivity));

            }
            return get_net_handle_from_values(get_record_offset(net_handle), connectivity, 
                                      get_handle_type(net_handle), get_node_record_offset(net_handle)+1);
        }
        //Otherwise, if we are going left from the first node or right from the last node, we pretend this is a 
        //snarl and use get_next_child of the offset, which will check for the length of this record properly
    }


    //Otherwise we need to skip through snarl and trivial snarl records
    //These will be arranged:
    //
    //trivial_snarl_record | trivial snarl length | snarl length | snarl record | snarl length | trivial snarl length | trivial snarl record ... 
    if (go_left && get_record_offset(net_handle) == get_first_node_offset() ) {
        //IF we are going left out the start of the chain
        if (get_start_id() != get_end_id()) {
            //If we are going left from the first node of a non-looping chain
            return net_handle;
        } else {
            //IF we are going left from the first node of a looping chain
            tuple<size_t, bool, size_t> last_child = get_last_child_offset();
            bool is_reversed_in_parent = std::get<1>(last_child) ? false : 
                TrivialSnarlRecord(get_first_node_offset(),records).get_is_reversed_in_parent(TrivialSnarlRecord(get_first_node_offset(),records).get_node_count()-1);
            
            return get_net_handle_from_values(std::get<0>(last_child), 
                                  (is_reversed_in_parent ? START_END : END_START),
                                  std::get<1>(last_child) ? SNARL_HANDLE : NODE_HANDLE,
                                  std::get<2>(last_child));
        }
    } else if (!go_left && get_record_offset(net_handle) == std::get<0>(get_last_child_offset())) {
        //IF we are going right out the end of the chain
        if (get_start_id() != get_end_id()) {
            //If we are going right from the last node of a non-looping chain
            return net_handle;
        } else {
            //IF we are going right from the last node of a looping chain, then
            //return the first node going in
            bool is_reversed_in_parent = TrivialSnarlRecord(get_first_node_offset(),records).get_is_reversed_in_parent(0);
            return get_net_handle_from_values(get_first_node_offset(), 
                                  (is_reversed_in_parent ? START_START : END_END),
                                  NODE_HANDLE,
                                  0);
        }
    }

    //At this point, we know that we need to leave a snarl or a trivial snarl and that we're in the middle of the chain
    //So the we need to skip to the start of the previous thing or to the end of this thing
    size_t next_pointer =  get_record_offset(net_handle) + 
                            (go_left ? -(*records)->at(get_record_offset(net_handle)-2)-2 : (*records)->at(get_record_offset(net_handle)-1)+2);

    if (SnarlDistanceIndex::is_nonroot_nonsimple_snarl(SnarlDistanceIndex::get_record_type((*records)->at(next_pointer)))) {
        //If the next thing is a snarl, then just return the snarl going in the direction we just moved in
        return get_net_handle_from_values(next_pointer, (go_left ? END_START : START_END), SNARL_HANDLE); 
    } else if (SnarlDistanceIndex::is_simple_snarl(SnarlDistanceIndex::get_record_type((*records)->at(next_pointer)))) {
        //If the next thing is a snarl, then just return the snarl going in the direction we just moved in
        return get_net_handle_from_values(next_pointer, (go_left ? END_START : START_END), SNARL_HANDLE, 1); 
    } else{

        //Otherwise, this is a node (trivial snarl of nodes) handle
        TrivialSnarlRecord next_record(next_pointer, records);

        //The node in the trivial snarl is the last one if we're going left, the first one if we're going right
        size_t node_offset = go_left ? next_record.get_node_count()-1 : 0;

        //Is the next node (the first thing in the snarl if we're going right, last thing if we're going left) reversed
        bool next_is_reversed_in_parent = next_record.get_is_reversed_in_parent( node_offset);

        //Get the direction of the next handle
        //Handle will point in whichever direction we moved in
        connectivity_t connectivity = go_left == next_is_reversed_in_parent ? START_END : END_START;

        return get_net_handle_from_values(next_pointer, connectivity, NODE_HANDLE, node_offset);
    }
}

bool SnarlDistanceIndex::ChainRecord::for_each_child(const std::function<bool(const net_handle_t&)>& iteratee) const {

    if (get_record_handle_type() == NODE_HANDLE) {
        //If this is a node pretending to be a chain, just do it for the node
        return iteratee(get_net_handle_from_values(record_offset, START_END, NODE_HANDLE));
    }


    size_t first_node_offset = get_first_node_offset();
    bool rev_in_parent = TrivialSnarlRecord(first_node_offset, records).get_is_reversed_in_parent(0);
    net_handle_t first_child = get_net_handle_from_values(first_node_offset, rev_in_parent ? END_START : START_END, NODE_HANDLE, 0);
    net_handle_t current_child = get_net_handle_from_values(first_node_offset, rev_in_parent ? END_START : START_END, NODE_HANDLE, 0);
    bool is_first = true;

    while (true) {
        if (!is_first && current_child == first_child){
            //Don't look at the first node a second time
            return true;
        }

        bool result = iteratee(current_child);
        if (result == false) {
            return false;
        }
        net_handle_t next_child = get_next_child(current_child, false);
        if (current_child == next_child || get_start_endpoint(get_connectivity(next_child)) == get_end_endpoint(get_connectivity(next_child))) {
            return true;
        } else {
            current_child = next_child;
        }
        is_first = false;
    }
    return true;
}


SnarlDistanceIndex::ChainRecordWriter::ChainRecordWriter (size_t pointer, record_t type, size_t node_count, bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector>* records){
    record_offset = pointer;
    records = records;
    SnarlTreeRecordWriter::record_offset = pointer;
    SnarlTreeRecordWriter::records = records;
    ChainRecord::record_offset = pointer;
    ChainRecord::records = records;
    (*records)->resize((*records)->size() + CHAIN_RECORD_SIZE);
    set_node_count(node_count);
    set_record_type(type);
}


void SnarlDistanceIndex::ChainRecordWriter::set_node_count(size_t node_count) {

    (*records)->at(record_offset + CHAIN_NODE_COUNT_OFFSET) = node_count;
}
size_t SnarlDistanceIndex::ChainRecord::get_depth() const {
    return (*records)->at(record_offset+CHAIN_DEPTH_OFFSET) ;
}
void SnarlDistanceIndex::ChainRecordWriter::set_depth(size_t depth) {

    (*records)->at(record_offset + CHAIN_DEPTH_OFFSET) = depth;
}

//The offset of the last child, if it is a snarl, and if it can loop
void SnarlDistanceIndex::ChainRecordWriter::set_last_child_offset(size_t offset, bool is_snarl, bool loopable) {

    (*records)->at(record_offset + CHAIN_LAST_CHILD_OFFSET) = ((offset << 2) | (is_snarl<<1)) | loopable;
}

void SnarlDistanceIndex::ChainRecordWriter::set_distance_left_start(size_t distance) {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + CHAIN_DISTANCE_LEFT_START_OFFSET) = (distance == std::numeric_limits<size_t>::max() ? 0 : distance + 1);
}
void SnarlDistanceIndex::ChainRecordWriter::set_distance_right_start(size_t distance) {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + CHAIN_DISTANCE_RIGHT_START_OFFSET) = (distance == std::numeric_limits<size_t>::max() ? 0 : distance + 1);
}
void SnarlDistanceIndex::ChainRecordWriter::set_distance_left_end(size_t distance) {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + CHAIN_DISTANCE_LEFT_END_OFFSET) = (distance == std::numeric_limits<size_t>::max() ? 0 : distance + 1);
}
void SnarlDistanceIndex::ChainRecordWriter::set_distance_right_end(size_t distance) {
    if (!has_distances(get_record_type())) {
        throw runtime_error("error: trying to access get distance in a distanceless index");
    }
    (*records)->at(record_offset + CHAIN_DISTANCE_RIGHT_END_OFFSET) = (distance == std::numeric_limits<size_t>::max() ? 0 : distance + 1);
}

//Add a snarl to the end of the chain and return a SnarlRecordWriter pointing to it
SnarlDistanceIndex::SnarlRecordWriter SnarlDistanceIndex::ChainRecordWriter::add_snarl(size_t snarl_size, record_t type, size_t hhl_size, size_t previous_child_offset) {

    size_t snarl_record_size = SnarlRecord::record_size(type, snarl_size, hhl_size);

    
    //The last child we saw was a trivial snarl and we need to record its size
    (*records)->at(previous_child_offset-1) = TrivialSnarlRecord(previous_child_offset, records).get_record_size();
    size_t current_size = (*records)->size();
    (*records)->resize(current_size+1);
    (*records)->at(current_size) = TrivialSnarlRecord(previous_child_offset, records).get_record_size();
    
    size_t start_i = (*records)->size();
    (*records)->resize(start_i+1);
    (*records)->at(start_i) = snarl_record_size;
    (*records)->reserve(start_i + snarl_record_size);
    SnarlRecordWriter snarl_record(snarl_size, records, type, hhl_size);
    snarl_record.set_parent_record_offset(get_offset());
    start_i = (*records)->size();
    (*records)->resize(start_i+1);
    (*records)->at(start_i) = snarl_record_size;
    return snarl_record;
}
//Add a simple snarl to the end of the chain and return a SnarlRecordWriter pointing to it
SnarlDistanceIndex::SimpleSnarlRecordWriter SnarlDistanceIndex::ChainRecordWriter::add_simple_snarl(size_t snarl_size, record_t type, size_t previous_child_offset) {

    size_t snarl_record_size = SIMPLE_SNARL_RECORD_SIZE + 2*snarl_size;

    
    //The last child we saw was a trivial snarl and we need to record its size
    (*records)->at(previous_child_offset-1) = TrivialSnarlRecord(previous_child_offset, records).get_record_size();
    size_t current_size = (*records)->size();
    (*records)->resize(current_size+1);
    (*records)->at(current_size) = TrivialSnarlRecord(previous_child_offset, records).get_record_size();
    
    //Add another place in snarl_tree_records for the size of this snarl
    size_t start_i = (*records)->size();
    (*records)->resize(start_i+1);
    (*records)->at(start_i) = snarl_record_size;
    (*records)->reserve(start_i + snarl_record_size);
    SimpleSnarlRecordWriter snarl_record(snarl_size, records, type);
    snarl_record.set_parent_record_offset(get_offset());
    start_i = (*records)->size();
    (*records)->resize(start_i+1);
    (*records)->at(start_i) = snarl_record_size;
    return snarl_record;
}


size_t SnarlDistanceIndex::ChainRecordWriter::add_node(nid_t node_id, size_t node_length, bool is_reversed_in_parent,
    size_t prefix_sum, size_t forward_loop, size_t reverse_loop, size_t component, size_t max_prefix_sum, size_t previous_child_offset, bool new_record,
    bool include_distances) {
    if ((SnarlDistanceIndex::is_trivial_snarl(SnarlDistanceIndex::get_record_type((*records)->at(previous_child_offset))))
            && (TrivialSnarlRecord(previous_child_offset, records).get_node_count() == MAX_TRIVIAL_SNARL_NODE_COUNT 
                || new_record || reverse_loop == 0)) {
        //If the last thing was a trivial snarl and it is full, then finish it off
        (*records)->at(previous_child_offset-1) = TrivialSnarlRecord(previous_child_offset, records).get_record_size();
        size_t current_size = (*records)->size();
        (*records)->resize(current_size+1);
        (*records)->at(current_size) = TrivialSnarlRecord(previous_child_offset, records).get_record_size();
    }

    if (previous_child_offset == 0
            || SnarlDistanceIndex::is_nonroot_nontrivial_snarl(SnarlDistanceIndex::get_record_type((*records)->at(previous_child_offset))) 
            || TrivialSnarlRecord(previous_child_offset, records).get_node_count() == MAX_TRIVIAL_SNARL_NODE_COUNT
            || reverse_loop == 0 || new_record) {
        //If the last thing was a snarl or nothing (previous_child_offset == 0, meaning that this is the 
        //first thing in the chain), then create a new trivial snarl
        size_t start_i = (*records)->size();
        (*records)->resize(start_i+1);
        (*records)->at(start_i) = 0; //Add a place for the size of this trivial snarl, to be filled in later

        //Set the node to point to the correct node record
        (*records)->at(get_node_pointer_offset(node_id,
                                              (*records)->at(MIN_NODE_ID_OFFSET),
                                              (*records)->at(COMPONENT_COUNT_OFFSET))) 
                    = (*records)->size();
        (*records)->at(get_node_pointer_offset(node_id,
                                              (*records)->at(MIN_NODE_ID_OFFSET),
                                              (*records)->at(COMPONENT_COUNT_OFFSET))+1) 
                    = 0;
    

        //Create a new trivial snarl record
        TrivialSnarlRecordWriter trivial_snarl_record((*records)->size(), include_distances ? DISTANCED_TRIVIAL_SNARL : TRIVIAL_SNARL, records, true);
        trivial_snarl_record.set_parent_record_offset(record_offset);
        trivial_snarl_record.set_node_count(1);
        if (include_distances) {
            trivial_snarl_record.set_prefix_sum(prefix_sum);
            trivial_snarl_record.set_max_prefix_sum(max_prefix_sum);;
            trivial_snarl_record.set_forward_loop(forward_loop);
            trivial_snarl_record.set_reverse_loop(reverse_loop);
            trivial_snarl_record.set_chain_component(component);
        }

        start_i = (*records)->size(); 
        if (include_distances) {
            (*records)->resize(start_i+2);
            (*records)->at(start_i) = (node_id<<1) | is_reversed_in_parent;
            (*records)->at(start_i+1) = node_length;
        } else {
            (*records)->resize(start_i+1);
            (*records)->at(start_i) = (node_id<<1) | is_reversed_in_parent;
        }

        //Return the offset offset of the new trivial snarl record
        return trivial_snarl_record.record_offset;
    } else {
        //The last thing on the chain is a trivial snarl that we are adding to the end of
        //
        //Add the node and its right prefix sum
        size_t start_i = (*records)->size(); 
        if (include_distances) {
            (*records)->resize(start_i+2);
            (*records)->at(start_i) = (node_id<<1) | is_reversed_in_parent;
            (*records)->at(start_i+1) = (*records)->at(start_i-1) + node_length;
        } else {
            (*records)->resize(start_i+1);
            (*records)->at(start_i) = (node_id<<1) | is_reversed_in_parent;
        }

        TrivialSnarlRecordWriter trivial_snarl_record(previous_child_offset, include_distances ? DISTANCED_TRIVIAL_SNARL : TRIVIAL_SNARL, records, false);

        //Increment the node count
        size_t old_node_count = trivial_snarl_record.get_node_count();
        trivial_snarl_record.set_node_count(old_node_count+1);
        //The forward loop value changes
        if (include_distances) {
            trivial_snarl_record.set_forward_loop(forward_loop);
        }

        //Set the node to point to the correct node record
        (*records)->at(get_node_pointer_offset(node_id,
                                              (*records)->at(MIN_NODE_ID_OFFSET),
                                              (*records)->at(COMPONENT_COUNT_OFFSET))) 
                    = previous_child_offset;
        (*records)->at(get_node_pointer_offset(node_id,
                                              (*records)->at(MIN_NODE_ID_OFFSET),
                                              (*records)->at(COMPONENT_COUNT_OFFSET))+1) 
                    = old_node_count ;
        //REturn the offset of the old trivial snarl record that we added t:
        return trivial_snarl_record.record_offset;
    }
}
string SnarlDistanceIndex::net_handle_as_string(const net_handle_t& net) const {
    net_handle_record_t type = get_handle_type(net);
    SnarlTreeRecord record (net, &snarl_tree_records);
    net_handle_record_t record_type = record.get_record_handle_type();
    string result = stringify(type) + " ";
    if (type == ROOT_HANDLE) {
        if (is_root_snarl(record.get_record_type())) {
            result += "root snarl";
            return result;
        } else {
            result += "root";
            return result;
        }
    } else if (type == NODE_HANDLE) {
        if (ends_at(net) == starts_at(net)) {
            result += "node" + std::to_string( node_id(net)) + (ends_at(net) == START ? "rev" : "fd") + " that is the end node of a looping chain";
            return result;
        }
        result += "node " + std::to_string( node_id(net)) + (ends_at(net) == START ? "rev" : "fd");
        return result;
    } else if (type == SNARL_HANDLE) {
        if (record.get_record_type() == ROOT) {
            result += "root snarl";
            return result;
        }
        if (get_node_record_offset(net) == 1) {
            result += "simple snarl ";
        } else if (get_node_record_offset(net) > 1) {
            result += " node in simple snarl ";
            result += node_id(net);
        } else {
            result += "snarl ";         
        }
    } else if (type == CHAIN_HANDLE && record_type == NODE_HANDLE) {
        result += "node " + std::to_string( NodeRecord(net, &snarl_tree_records).get_node_id())
               + (ends_at(net) == START ? "rev" : "fd") + " pretending to be a chain";
        return result;
    }  else if (type == CHAIN_HANDLE && record_type == SNARL_HANDLE) {
        result += "node " + std::to_string( SimpleSnarlRecord(net, &snarl_tree_records).get_node_id())
               + (ends_at(net) == START ? "rev" : "fd") + " pretending to be a chain in a simple snarl";
        return result;
    }else if (type == CHAIN_HANDLE) {
        result += "chain ";
    } else if (type == SENTINEL_HANDLE) {
        result += "sentinel of snarl ";
    } else {
        throw runtime_error("error: Unknown net_handle_t type");
    }
    result += ( std::to_string(record.get_start_id())
            + (record.get_start_orientation() ? "rev" : "fd")
            + "->"
            + std::to_string(record.get_end_id())
            + (record.get_end_orientation() ? "rev" : "fd"));
    result += " traversing ";
    result += (starts_at(net) == START ? "start" : (starts_at(net) == END ? "end" : "tip"));
    result += "->";
    result += (ends_at(net) == START ? "start" : (ends_at(net) == END ? "end" : "tip"));
    return result;
}

bool SnarlDistanceIndex::traverse_decomposition(const std::function<bool(const net_handle_t&)>& snarl_iteratee,
                                                const std::function<bool(const net_handle_t&)>& chain_iteratee,
                                                const std::function<bool(const net_handle_t&)>& node_iteratee) const {
    net_handle_t root = get_root();
    return for_each_child(root, [&] (const net_handle_t& child) {
        return traverse_decomposition_helper(child, snarl_iteratee, chain_iteratee, node_iteratee);
    });

}
bool SnarlDistanceIndex::traverse_decomposition_helper(const net_handle_t& net,
                                                const std::function<bool(const net_handle_t&)>& snarl_iteratee,
                                                const std::function<bool(const net_handle_t&)>& chain_iteratee,
                                                const std::function<bool(const net_handle_t&)>& node_iteratee) const {
    if (is_node(net)) {
        //If this is a node, then run iteratee and we're done
        return node_iteratee(net);
    } else if (is_snarl(net)) {
        //IF this is a snarl, run iteratee and if it returns false, return false to stop iterating
        bool result = snarl_iteratee(net);
        if (!result){
            return false;
        }
    } else if (is_chain(net)) {
        //IF this is a chain, run iteratee and if it returns false, return false to stop iterating
        bool result = chain_iteratee(net);
        if (!result){
            return false;
        }
    } else {
        throw runtime_error("error: traversing the wrong kind of net_handle_t in traverse_decomposition");
    }
    //If this was a snarl or chain, recurse
    return for_each_child(net, [&] (const net_handle_t& child) {
        return traverse_decomposition_helper(child, snarl_iteratee, chain_iteratee, node_iteratee);
    });

}



//Print entire index to cout as csv:
//self, parent, # children, depth
void SnarlDistanceIndex::print_self() const {
    cout << "#self,parent,#children,depth" << endl; 
    print_descendants_of(get_root());

}
void SnarlDistanceIndex::print_descendants_of(const net_handle_t net) const {
    SnarlTreeRecord record (net, &snarl_tree_records);
    //What the record thinks it is
    net_handle_record_t record_type = record.get_record_handle_type();
    if (record_type == NODE_HANDLE) {
        net_handle_t parent = get_parent(net);
        cout << net_handle_as_string(net) << "," << net_handle_as_string(parent) << ",1,"<< get_depth(net) << endl; 
        return;
    } else {
        size_t child_count;
        string parent;
        if (record_type == ROOT_HANDLE) {
            parent = "none";
            child_count = RootRecord(get_root(), &snarl_tree_records).get_connected_component_count();
        } else {
            parent = net_handle_as_string(get_parent(net));
            if (record_type == CHAIN_HANDLE) {
                child_count = ChainRecord(net, &snarl_tree_records).get_node_count();
            } else if (is_nonroot_nonsimple_snarl(record.get_record_type())) {
 
                child_count = SnarlRecord(net, &snarl_tree_records).get_node_count();
            } else if (is_trivial_snarl(record.get_record_type())) {
                child_count = TrivialSnarlRecord(get_record_offset(net), &snarl_tree_records).get_node_count();
            } else if (is_simple_snarl(record.get_record_type())) {
                child_count = SimpleSnarlRecord(net, &snarl_tree_records).get_node_count();
            } else {
                throw runtime_error("error: printing the wrong kind of record");
            }
                                
        }
        cout << net_handle_as_string(net) << "," << parent << "," << child_count << ","<< get_depth(net) << endl; 
        if (!is_node(net)) {
        for_each_child(net, [&](const net_handle_t& child) {
            print_descendants_of(child);
        });
        }
    }
}

void SnarlDistanceIndex::print_snarl_stats() const {
    cout << "#start_id\tend_id\tsnarl_size\tsnarl_depth" << endl;
    traverse_decomposition(
        [&](const net_handle_t& snarl_child) {
            //Iteratee for a snarl child
            SnarlTreeRecord record(snarl_child, &snarl_tree_records);

            //Get the number of children depending on the type of record
            size_t child_count;
            if (is_nonroot_nonsimple_snarl(record.get_record_type())) {

                child_count = SnarlRecord(snarl_child, &snarl_tree_records).get_node_count();
            } else if (is_simple_snarl(record.get_record_type())) {
                child_count = SimpleSnarlRecord(snarl_child, &snarl_tree_records).get_node_count();
            } else {
                throw runtime_error("error: getting the snarl child count of the wrong type of record");
            }
            //Print the stats
            cout << record.get_start_id() << "\t" 
                 << record.get_end_id() << "\t" 
                 << child_count << "\t" 
                 << get_depth(snarl_child) << endl; 
            return true;
        }, 
        [&](const net_handle_t& chain_child) {
            //Iteratee of a chain child -do nothing
            return true;
        }, 
        [&](const net_handle_t& node_child) {
            //Iteratee of a node child - do nothing
            return true;
        });
}

void SnarlDistanceIndex::write_snarls_to_json() const {
    traverse_decomposition(
        [&](const net_handle_t& snarl_child) {
            //Iteratee for a snarl child

            //Make a new json object for this snarl
            json_t* out_json = json_object();
            json_object_set_new(out_json, "type", json_string("snarl")); 
            SnarlTreeRecord record(snarl_child, &snarl_tree_records);

            //Get the start node pointing in
            net_handle_t start_bound = get_node_from_sentinel(get_bound(snarl_child, false, true));
            json_t* start_json = json_object();
            json_object_set_new(start_json, "node_id", json_integer(node_id(start_bound)));
            json_object_set_new(start_json, "is_rev", json_boolean(ends_at(start_bound) == START));
            //Set the start node
            json_object_set_new(out_json, "start_node", start_json);

            //Get the end node pointing in
            net_handle_t end_bound = get_node_from_sentinel(get_bound(snarl_child, true, true));
            json_t* end_json = json_object();
            json_object_set_new(end_json, "node_id", json_integer(node_id(end_bound)));
            json_object_set_new(end_json, "is_rev", json_boolean(ends_at(end_bound) == START));
            //Set the end node
            json_object_set_new(out_json, "end_node", end_json);

            //Set the parent
            net_handle_t parent = get_parent(snarl_child);
                json_t* parent_json = json_object();
            if (is_root(parent)) {
                json_object_set_new(parent_json, "type", json_string("root"));
            } else {
                json_object_set_new(parent_json, "type", json_string("chain"));
                net_handle_t parent_start = get_bound(parent, false, true);
                json_object_set_new(parent_json, "start_id", json_integer(node_id(parent_start)));
                net_handle_t parent_end = get_bound(parent, true, false);
                json_object_set_new(parent_json, "end_id", json_integer(node_id(parent_end)));
            }
            json_object_set_new(out_json, "parent", parent_json);

            //Get the number of children depending on the type of record
            if (is_nonroot_nonsimple_snarl(record.get_record_type())) {

                size_t child_count = SnarlRecord(snarl_child, &snarl_tree_records).get_node_count();
                json_object_set_new(out_json, "child_count", json_integer(child_count)); 
            } else if (is_simple_snarl(record.get_record_type())) {
                size_t child_count = SimpleSnarlRecord(snarl_child, &snarl_tree_records).get_node_count();
                json_object_set_new(out_json, "child_count", json_integer(child_count)); 
            } else {
                throw runtime_error("error: getting the snarl child count of the wrong type of record");
            }

            //Set the min and max length, if applicable
            if (has_distances(record.get_record_type()) && is_nonroot_nontrivial_snarl(record.get_record_type()) && !is_oversized_snarl(record.get_record_type())){
                // TODO: May need a new accessor to match just snarls and simple snarls that have distances.
                json_object_set_new(out_json, "minimum_length", json_integer(minimum_length(snarl_child))); 
                json_object_set_new(out_json, "maximum_length", json_integer(maximum_length(snarl_child))); 
            }

            //Set the depth
            json_object_set_new(out_json, "depth", json_integer(get_depth(snarl_child)));

            cout << json_dumps(out_json, JSON_ENCODE_ANY) << endl;

            return true;
        }, 
        [&](const net_handle_t& chain_child) {
            //Iteratee of a chain child
            if (is_trivial_chain(chain_child)) {
                //If this is a trivial chain, do nothing, handle it as a node
                return true;
            }
            json_t* out_json = json_object();
            json_object_set_new(out_json, "type", json_string("chain"));

            //Make a new json object for this chain
            ChainRecord record(chain_child, &snarl_tree_records);

            //Get the start node pointing in
            net_handle_t start_bound = get_bound(chain_child, false, true);
            json_t* start_json = json_object();
            json_object_set_new(start_json, "node_id", json_integer(node_id(start_bound)));
            json_object_set_new(start_json, "is_rev", json_boolean(ends_at(start_bound) == START));
            //Set the start node
            json_object_set_new(out_json, "start_node", start_json);

            //Get the end node pointing in
            net_handle_t end_bound = get_bound(chain_child, true, true);
            json_t* end_json = json_object();
            json_object_set_new(end_json, "node_id", json_integer(node_id(end_bound)));
            json_object_set_new(end_json, "is_rev", json_boolean(ends_at(end_bound) == START));
            //Set the end node
            json_object_set_new(out_json, "end_node", end_json);

            //Set the parent
            net_handle_t parent = get_parent(chain_child);
                json_t* parent_json = json_object();
            if (is_root(parent)) {
                json_object_set_new(parent_json, "type", json_string("root"));
            } else {
                json_object_set_new(parent_json, "type", json_string("snarl"));
                net_handle_t parent_start = get_bound(parent, false, true);
                json_object_set_new(parent_json, "start_id", json_integer(node_id(parent_start)));
                net_handle_t parent_end = get_bound(parent, true, false);
                json_object_set_new(parent_json, "end_id", json_integer(node_id(parent_end)));
            }
            json_object_set_new(out_json, "parent", parent_json);

            //Get the number of children depending on the type of record
            json_object_set_new(out_json, "child_count", json_integer(record.get_node_count())); 

            //Set the min and max length, if applicable
            if (record.get_record_type() == DISTANCED_CHAIN) {
                json_object_set_new(out_json, "minimum_length", json_integer(minimum_length(chain_child))); 
                json_object_set_new(out_json, "maximum_length", json_integer(maximum_length(chain_child))); 
            }

            //Set the depth
            json_object_set_new(out_json, "depth", json_integer(get_depth(chain_child)));

            cout << json_dumps(out_json, JSON_ENCODE_ANY) << endl;
            return true;
        }, 
        [&](const net_handle_t& node_child) {
            //Iteratee of a node child
            json_t* out_json = json_object();

            net_handle_t parent = get_parent(node_child);
            if (is_trivial_chain(parent)) {
                json_object_set_new(out_json, "type", json_string("trivial_chain"));
                json_object_set_new(out_json, "node_id", json_integer(node_id(node_child)));
                net_handle_t grandparent_snarl = get_parent(parent);

                json_t* parent_json = json_object();
                if (is_root(parent)) {
                    json_object_set_new(parent_json, "type", json_string("root"));
                } else {
                    json_object_set_new(parent_json, "type", json_string("snarl"));
                    net_handle_t parent_start = get_bound(grandparent_snarl, false, true);
                    json_object_set_new(parent_json, "start_id", json_integer(node_id(parent_start)));
                    net_handle_t parent_end = get_bound(grandparent_snarl, true, false);
                    json_object_set_new(parent_json, "end_id", json_integer(node_id(parent_end)));
                }
                json_object_set_new(out_json, "parent", parent_json);
            } else {
                //Otherwise this is just a node
                json_object_set_new(out_json, "type", json_string("node"));
                json_object_set_new(out_json, "node_id", json_integer(node_id(node_child)));

                json_t* parent_json = json_object();
                if (is_root(parent)) {
                    json_object_set_new(parent_json, "type", json_string("root"));
                } else {
                    json_object_set_new(parent_json, "type", json_string("chain"));
                    net_handle_t parent_start = get_bound(parent, false, true);
                    json_object_set_new(parent_json, "start_id", json_integer(node_id(parent_start)));
                    net_handle_t parent_end = get_bound(parent, true, false);
                    json_object_set_new(parent_json, "end_id", json_integer(node_id(parent_end)));
                }
                json_object_set_new(out_json, "parent", parent_json);
            }
            json_object_set_new(out_json, "length", json_integer(minimum_length(node_child)));
            cout << json_dumps(out_json, JSON_ENCODE_ANY) << endl;
            return true;
        });
}

std::tuple<size_t, size_t, size_t> SnarlDistanceIndex::get_usage() {
    return snarl_tree_records.get_usage();
}


void SnarlDistanceIndex::get_snarl_tree_records(const vector<TemporaryDistanceIndex*>& temporary_indexes, const HandleGraph* graph) {


    //Convert temporary distance indexes into the final index stored as a single vector
    size_t total_component_count = 0;
    handlegraph::nid_t min_node_id = 0;
    handlegraph::nid_t max_node_id = 0;

    size_t maximum_distance = 0;
    size_t maximum_index_size = 0;
    size_t maximum_tree_depth = 0;
    // Start with 10
    //TODO: This should maybe be a const static variable somewhere
    size_t max_bits = 10;
    //The maximum distance value that gets stored in the index, and the maximum possible
    //length of the index

    /*Go through each of the indexes to count how many nodes, components, etc */
    for (const TemporaryDistanceIndex* temp_index : temporary_indexes) {
        total_component_count += temp_index->root_structure_count;
        min_node_id = min_node_id == 0 ? temp_index->min_node_id
                                       : std::min(min_node_id, temp_index->min_node_id);
        max_node_id = std::max(max_node_id, temp_index->max_node_id);
        maximum_distance = std::max(temp_index->max_distance, maximum_distance);
        maximum_index_size += temp_index->get_max_record_length();
        maximum_tree_depth = std::max(maximum_tree_depth, temp_index->max_tree_depth);
        max_bits = std::max(max_bits, temp_index->max_bits);
    }



    //Set the width of the values to the minimum needed to fit everything
    size_t max_dist_bit_width = bit_width(std::max(maximum_distance, (size_t) max_node_id));
    size_t max_address_bit_width = bit_width(maximum_index_size);
    size_t new_width = std::max(std::max(max_dist_bit_width, max_address_bit_width)+2, max_bits);
    if (new_width > 64) {
        cerr << "The bit width for the distance index is greater than 64 bits, which may cause problems" << endl;
    }
    snarl_tree_records->width(new_width);

    snarl_tree_records->reserve(maximum_index_size);

    /*Allocate memory for the root and the nodes */
    RootRecordWriter root_record(0, total_component_count, max_node_id-min_node_id+1, maximum_tree_depth, min_node_id, &snarl_tree_records);

    /*Now go through each of the chain/snarl indexes and copy them into snarl_tree_records
     * Walk down the snarl tree and fill in children
     */
    // maps <index into temporary_indexes, <record type, index into chain/snarl/node records>> to new offset
    unordered_map<std::pair<size_t, temp_record_ref_t>, size_t> record_to_offset;
    //Set the root index
    for (size_t temp_index_i = 0 ; temp_index_i < temporary_indexes.size() ; temp_index_i++) {
        //Any root will point to the same root
        record_to_offset.emplace(make_pair(temp_index_i,make_pair(TEMP_ROOT, 0)), 0);
    }
    //Go through each separate temporary index, corresponding to separate connected components
    //Root components of all the indexes are numbered consecutively
    size_t first_component_num = 0;
    for (size_t temp_index_i = 0 ; temp_index_i < temporary_indexes.size() ; temp_index_i++) {
        TemporaryDistanceIndex* temp_index = temporary_indexes[temp_index_i];
        //Get a stack of temporary snarl tree records to be added to the index
        //Initially, it contains only the root components
        //This reverses the order of the connected components but I don't think that matters
        vector<temp_record_ref_t> temp_record_stack = temp_index->components;

        while (!temp_record_stack.empty()) {
            temp_record_ref_t current_record_index = temp_record_stack.back();
            temp_record_stack.pop_back();


            if (current_record_index.first == TEMP_CHAIN) {
                /*Add a new chain to the index. Each of the chain's child snarls and nodes will also
                 * be added here
                 */
                const TemporaryDistanceIndex::TemporaryChainRecord& temp_chain_record = temp_index->get_chain(current_record_index);
                if (!temp_chain_record.is_trivial) {
                    //If this chain contains at least two nodes
                    record_to_offset.emplace(make_pair(temp_index_i,current_record_index), snarl_tree_records->size());


                    //If this chain's parent is a root snarl, then the temporary parent is a snarl but we consider 
                    //the parent to be the root and the chain to have depth 1
                    bool is_child_of_root_snarl = false;
                    if (temp_chain_record.parent.first == TEMP_SNARL) {
                        const TemporaryDistanceIndex::TemporarySnarlRecord& temp_parent_record = temp_index->get_snarl(temp_chain_record.parent);
                        if (temp_parent_record.is_root_snarl) {
                            is_child_of_root_snarl = true;
                        }
                    }

                    //We don't keep track of distances for this chain if the snarl_size_limit is 0, or if we only want top-level chain distances
                    // and this is a nested chain
                    bool ignore_distances = (snarl_size_limit == 0) || 
                                            (only_top_level_chain_distances && !(temp_chain_record.parent.first == TEMP_ROOT || is_child_of_root_snarl));

                    ChainRecordWriter chain_record_constructor;

                    if (temp_chain_record.chain_components.back() == 0 || ignore_distances) {
                        record_t record_type = encode_chain(!ignore_distances, false);
                        chain_record_constructor = ChainRecordWriter(snarl_tree_records->size(), record_type,
                                                               temp_chain_record.prefix_sum.size(), &snarl_tree_records);
                        chain_record_constructor.set_start_end_connected();
                    } else {
                        chain_record_constructor = ChainRecordWriter(snarl_tree_records->size(), encode_chain(!ignore_distances, true),
                                                               temp_chain_record.prefix_sum.size(), &snarl_tree_records);
                    }
                    chain_record_constructor.set_parent_record_offset(
                            record_to_offset[make_pair(temp_index_i, temp_chain_record.parent)]);
                    if (!ignore_distances) {
                        chain_record_constructor.set_distance_left_start(temp_chain_record.distance_left_start);
                        chain_record_constructor.set_distance_right_start(temp_chain_record.distance_right_start);
                        chain_record_constructor.set_distance_left_end(temp_chain_record.distance_left_end);
                        chain_record_constructor.set_distance_right_end(temp_chain_record.distance_right_end);
                    }

                    //Set the depth of this chain
                    if (temp_chain_record.parent.first == TEMP_ROOT || is_child_of_root_snarl) {
                        //If the parent is the root, then its depth is 0
                        chain_record_constructor.set_depth(1);
                    } else {
                        //Otherwise, go to the grandparent chain and add 1
                        size_t parent_record_offset = chain_record_constructor.get_parent_record_offset();
                        size_t grandparent_record_offset = SnarlRecord(parent_record_offset, &snarl_tree_records).get_parent_record_offset();  
                        chain_record_constructor.set_depth(
                            ChainRecord(grandparent_record_offset, &snarl_tree_records).get_depth() + 1);
                    }


                    if (!ignore_distances) {
                        chain_record_constructor.set_min_length(temp_chain_record.min_length);
                        chain_record_constructor.set_max_length(temp_chain_record.max_length);
                    }
                    chain_record_constructor.set_rank_in_parent(temp_chain_record.rank_in_parent);
                    chain_record_constructor.set_start_node(temp_chain_record.start_node_id, temp_chain_record.start_node_rev);
                    chain_record_constructor.set_end_node(temp_chain_record.end_node_id, temp_chain_record.end_node_rev);


                    size_t chain_node_i = 0; //How far along the chain are we?
                    pair<size_t, bool> last_child_offset = make_pair(0, false);

                    //A trivial snarl might not actually be a trivial snarl if it has loops in it. Keep track of whether the last 
                    //child node actually had a 0 value forward loop, making the last child a non-trivial snarl
                    bool last_child_was_nontrivial_snarl = false;

                    for (size_t child_record_index_i = 0 ; child_record_index_i < temp_chain_record.children.size() ; child_record_index_i++) {
                        const temp_record_ref_t& child_record_index = temp_chain_record.children[child_record_index_i];
                        //Go through each node and snarl in the chain and add them to the index

                        if (child_record_index.first == TEMP_NODE) {
                            //Add a node to the chain
                            if (!(chain_node_i != 0 && child_record_index == temp_chain_record.children.front())) {
                                //If this is not a looping chain, then we haven't seen this node yet

                                //Get the temporary node record
                                const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = temp_index->get_node(child_record_index); 

                                //Make a new node record
                                size_t new_offset = chain_record_constructor.add_node(
                                        temp_node_record.node_id, temp_node_record.node_length, temp_node_record.reversed_in_parent,
                                        temp_chain_record.prefix_sum[chain_node_i], temp_chain_record.forward_loops[chain_node_i],
                                        temp_chain_record.backward_loops[chain_node_i], temp_chain_record.chain_components[chain_node_i],
                                        temp_chain_record.max_prefix_sum[chain_node_i],
                                        last_child_offset.first, last_child_was_nontrivial_snarl,
                                        !ignore_distances);

                                //Remember this node as the last thing in the chain
                                last_child_offset = make_pair(new_offset, false);
                                if (temp_chain_record.forward_loops[chain_node_i] == 0) {
                                    last_child_was_nontrivial_snarl = true;
                                } else {
                                    last_child_was_nontrivial_snarl = false;
                                }

                            } 

                            chain_node_i++;


                        } else {
                            //Add a snarl to the chain
                            assert(child_record_index.first == TEMP_SNARL);
                            //Get the temporary snarl record
                            TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index->get_snarl(child_record_index);
                            if (!temp_snarl_record.is_trivial && !temp_snarl_record.is_simple) {
                                //If this is an actual snarl that we need to make


                                //Add the snarl to the chain, and get back the record to fill it in

                                bool ignore_distances = (snarl_size_limit == 0) || only_top_level_chain_distances;
                                
                                record_t record_type = encode_nonroot_nonsimple_snarl(!ignore_distances, temp_snarl_record.is_regular, !ignore_distances && temp_snarl_record.node_count > snarl_size_limit);

                                SnarlRecordWriter snarl_record_constructor =
                                    chain_record_constructor.add_snarl(temp_snarl_record.node_count, record_type, temp_snarl_record.hub_labels.size(), last_child_offset.first);

                                //Record how to find the new snarl record
                                record_to_offset.emplace(make_pair(temp_index_i, child_record_index), snarl_record_constructor.record_offset);

                                //Fill in snarl info
                                if (!ignore_distances) {
                                    snarl_record_constructor.set_min_length(temp_snarl_record.min_length);
                                    snarl_record_constructor.set_max_length(temp_snarl_record.max_length);
                                }
                                snarl_record_constructor.set_distance_start_start(temp_snarl_record.distance_start_start);
                                snarl_record_constructor.set_distance_end_end(temp_snarl_record.distance_end_end);

                                //Add distances and record connectivity
                                if (is_oversized_snarl(record_type)) {
                                    // We need to copy the packed hub label vector into place.
                                    for (size_t i = 0; i < temp_snarl_record.hub_labels.size(); i++) {
                                        // TODO: Make this an std::copy or something.
                                        snarl_record_constructor.set_hhl_entry(i, temp_snarl_record.hub_labels.at(i));
                                    }
                                    // TODO: When should we call
                                    // snarl_record_constructor.set_tip_tip_connected()?
                                    // Add code to determine that somewhere!
                                } else {
                                    // Store individual distance entries.
                                    if (!ignore_distances) {
                                        temp_snarl_record.for_each_distance([&](size_t rank1, bool right_side1, size_t rank2, bool right_side2, size_t distance) {
                                            snarl_record_constructor.set_distance(rank1, right_side1, rank2, right_side2, distance);

                                            if (temp_snarl_record.tippy_child_ranks.count(rank1)
                                                && temp_snarl_record.tippy_child_ranks.count(rank2)) {
                                                snarl_record_constructor.set_tip_tip_connected();
                                            }
                                        });
                                    }
                                }
                                // The final record has the distances now, so free the temporary copies
                                vector<size_t>().swap(temp_snarl_record.distances);
                                vector<size_t>().swap(temp_snarl_record.hub_labels);
                                //Now set the connectivity of this snarl
                                if (temp_snarl_record.distance_start_start != std::numeric_limits<size_t>::max()) {
                                    snarl_record_constructor.set_start_start_connected();
                                }
                                if (temp_snarl_record.min_length != std::numeric_limits<size_t>::max()) {
                                    snarl_record_constructor.set_start_end_connected();
                                }
                                if (temp_snarl_record.distance_end_end != std::numeric_limits<size_t>::max()) {
                                    snarl_record_constructor.set_end_end_connected();
                                }

                                for (const temp_record_ref_t& child : temp_snarl_record.children) {
                                    temp_record_stack.emplace_back(child);
                                }


                                //Add connectivity in chain based on this snarl
                                //TODO: Tip-tip connected?
                                //TODO: What about if it's a multicomponent chain?
                                if (!snarl_record_constructor.get_is_reversed_in_parent()) {
                                    //If this snarl is oriented forward in the chain
                                    if (snarl_record_constructor.is_start_start_connected()) {
                                        chain_record_constructor.set_start_start_connected();
                                    }
                                    if (snarl_record_constructor.is_end_end_connected()) {
                                        chain_record_constructor.set_end_end_connected();
                                    }
                                    if (snarl_record_constructor.is_start_tip_connected()) {
                                        chain_record_constructor.set_start_tip_connected();
                                    }
                                    if (snarl_record_constructor.is_end_tip_connected()) {
                                        chain_record_constructor.set_end_tip_connected();
                                    }
                                } else {
                                    //If this snarl is oriented backwards in the chain
                                    if (snarl_record_constructor.is_start_start_connected()) {
                                        chain_record_constructor.set_end_end_connected();
                                    }
                                    if (snarl_record_constructor.is_end_end_connected()) {
                                        chain_record_constructor.set_start_start_connected();
                                    }
                                    if (snarl_record_constructor.is_start_tip_connected()) {
                                        chain_record_constructor.set_end_tip_connected();
                                    }
                                    if (snarl_record_constructor.is_end_tip_connected()) {
                                        chain_record_constructor.set_start_tip_connected();
                                    }
                                }
                                last_child_offset = make_pair(snarl_record_constructor.record_offset, true);
                            } else if (!temp_snarl_record.is_trivial && temp_snarl_record.is_simple) {
                                //Make a simple snarl

                                //Add the snarl to the chain, and get back the record to fill it in
                                bool ignore_distances = (snarl_size_limit == 0) || only_top_level_chain_distances;

                                record_t record_type = encode_simple_snarl(!ignore_distances);

                                SimpleSnarlRecordWriter snarl_record_constructor =
                                    chain_record_constructor.add_simple_snarl(temp_snarl_record.node_count, record_type, last_child_offset.first);

                                //Record how to find the new snarl record
                                record_to_offset.emplace(make_pair(temp_index_i, child_record_index), snarl_record_constructor.record_offset);

                                //Fill in snarl info
                                if (!ignore_distances) {
                                    snarl_record_constructor.set_min_length(temp_snarl_record.min_length);
                                    snarl_record_constructor.set_max_length(temp_snarl_record.max_length);
                                }

                                //Add the children of the simple snarl
                                for (size_t i = 0 ; i < temp_snarl_record.node_count ; i++ ) {
                                    const temp_record_ref_t& child_index = temp_snarl_record.children[i];
                                    if( child_index.first == TEMP_CHAIN) {
                                        if (temp_index->get_chain(child_index).children.size() != 1) {
                                           throw runtime_error("size of children != 1, child index: "+to_string(child_index.second)+", bounding nodes: "+to_string(temp_snarl_record.start_node_id)+" "+to_string(temp_snarl_record.end_node_id));
                                        }
                                        const temp_record_ref_t& node_index = temp_index->get_chain(child_index).children.front();
                                        const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = temp_index->get_node(node_index);
                                        //If there is a way to go from the node forward to the start node,
                                        //then it is reversed
                                        size_t rank = temp_index->get_chain(child_index).rank_in_parent;

                                        snarl_record_constructor.add_child(i+2, temp_node_record.node_id,  
                                                temp_node_record.node_length, temp_node_record.reversed_in_parent);
                                    } else {
                                        assert(child_index.first == TEMP_NODE);
                                        const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = temp_index->get_node(child_index);
                                        size_t rank =temp_node_record.rank_in_parent;
                                        snarl_record_constructor.add_child(i+2, temp_node_record.node_id,  
                                                temp_node_record.node_length, temp_node_record.reversed_in_parent);
                                    }
                                }

                                last_child_offset = make_pair(snarl_record_constructor.record_offset, true);
                            }
                            last_child_was_nontrivial_snarl = false;
                        }
                    }
                    //Does the chain loop and is the last node connected to the rest of the chain through the last snarl
                    bool last_node_connected = temp_chain_record.loopable && (temp_chain_record.start_node_id==temp_chain_record.end_node_id);
                    chain_record_constructor.set_last_child_offset(last_child_offset.first, last_child_offset.second, last_node_connected);
                    //Finish the chain by adding two 0's
                } else {
                    //If the chain is trivial, then only record the node
                    assert(temp_chain_record.children.size() == 1);
                    assert(temp_chain_record.children[0].first == TEMP_NODE);
                    const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = temp_index->get_node(temp_chain_record.children[0]);


                    bool ignore_distances = (snarl_size_limit == 0) || only_top_level_chain_distances;

                    record_t record_type = encode_node(!ignore_distances);
                    NodeRecordWriter node_record(snarl_tree_records->size(), 0, record_type, &snarl_tree_records, temp_node_record.node_id);
                    node_record.set_node_id(temp_node_record.node_id);
                    node_record.set_rank_in_parent(temp_chain_record.rank_in_parent);
                    node_record.set_parent_record_offset(record_to_offset[make_pair(temp_index_i, temp_chain_record.parent)]);
                    if (!ignore_distances) {
                        node_record.set_node_length(temp_node_record.node_length);
                        node_record.set_distance_left_start(temp_chain_record.distance_left_start);
                        node_record.set_distance_right_start(temp_chain_record.distance_right_start);
                        node_record.set_distance_left_end(temp_chain_record.distance_left_end);
                        node_record.set_distance_right_end(temp_chain_record.distance_right_end);
                    }

                    record_to_offset.emplace(make_pair(temp_index_i, current_record_index), node_record.record_offset);

                }
            } else if (current_record_index.first == TEMP_SNARL) {

                bool ignore_distances = (snarl_size_limit == 0) || only_top_level_chain_distances;
                TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index->get_snarl(current_record_index);
                bool is_oversized = !temp_snarl_record.hub_labels.empty();
                //This is a root-level snarl
                record_t record_type = encode_root_snarl(!ignore_distances, is_oversized);

                record_to_offset.emplace(make_pair(temp_index_i,current_record_index), snarl_tree_records->size());

                SnarlRecordWriter snarl_record_constructor (temp_snarl_record.node_count, &snarl_tree_records, record_type, temp_snarl_record.hub_labels.size());

                //Fill in snarl info
                snarl_record_constructor.set_parent_record_offset(0);

                //Add distances and record connectivity

                if (!ignore_distances) {
                    if (is_oversized) {
                        //Copy the packed hub-label vector into the record.
                        for (size_t i = 0; i < temp_snarl_record.hub_labels.size(); i++) {
                            snarl_record_constructor.set_hhl_entry(i, temp_snarl_record.hub_labels.at(i));
                        }
                    } else {
                        temp_snarl_record.for_each_distance([&](size_t rank1, bool right_side1, size_t rank2, bool right_side2, size_t distance) {
                            snarl_record_constructor.set_distance(rank1, right_side1, rank2, right_side2, distance);
                        });
                    }
                }
                // The final record has the distances now, so free the temporary copies
                vector<size_t>().swap(temp_snarl_record.distances);
                vector<size_t>().swap(temp_snarl_record.hub_labels);

                for (const temp_record_ref_t& child : temp_snarl_record.children) {
                        temp_record_stack.emplace_back(child);
                }

            } else {
                assert(current_record_index.first == TEMP_NODE);
                //and then add them all after adding the snarl
                const TemporaryDistanceIndex::TemporaryNodeRecord& temp_node_record = temp_index->get_node(current_record_index);

                bool ignore_distances = (snarl_size_limit == 0) || only_top_level_chain_distances;
                record_t record_type = encode_node(!ignore_distances);
                NodeRecordWriter node_record(snarl_tree_records->size(), 0, record_type, &snarl_tree_records, temp_node_record.node_id);
                node_record.set_node_id(temp_node_record.node_id);
                node_record.set_rank_in_parent(temp_node_record.rank_in_parent);
                node_record.set_parent_record_offset(record_to_offset[make_pair(temp_index_i, temp_node_record.parent)]);

                if (!ignore_distances) {
                    node_record.set_node_length(temp_node_record.node_length);
                    node_record.set_distance_left_start(temp_node_record.distance_left_start);
                    node_record.set_distance_right_start(temp_node_record.distance_right_start);
                    node_record.set_distance_left_end(temp_node_record.distance_left_end);
                    node_record.set_distance_right_end(temp_node_record.distance_right_end);
                }

                record_to_offset.emplace(make_pair(temp_index_i, current_record_index), node_record.record_offset);
            }
        }

        for (size_t component_num = 0 ; component_num < temp_index->components.size() ; component_num++){
            const temp_record_ref_t& component_index = temp_index->components[component_num];
            //Let the root record know that it has another root
            root_record.add_component(first_component_num + component_num,record_to_offset[make_pair(temp_index_i,component_index)]);

            SnarlTreeRecord record (record_to_offset[make_pair(temp_index_i, component_index)],
                                    &snarl_tree_records);
            SnarlTreeRecordWriter record_constructor(record_to_offset[make_pair(temp_index_i, component_index)],
                                                              &snarl_tree_records);

            if (record.get_record_handle_type() != SNARL_HANDLE) {
                //Nodes, chains, and root snarls at the root get their global component number
                record_constructor.set_rank_in_parent(first_component_num + component_num);
            }
            if (!is_root_snarl(record.get_record_type())) {
                //If this isn't a root snarl
                handle_t start_out = graph->get_handle(record.get_start_id(), !record.get_start_orientation());
                handle_t end_out = graph->get_handle(record.get_end_id(), record.get_end_orientation());
                handle_t start_in = graph->get_handle(record.get_start_id(), record.get_start_orientation());
                handle_t end_in = graph->get_handle(record.get_end_id(), !record.get_end_orientation());


                graph->follow_edges(start_out, false, [&](const handle_t& h) {
                    if (h == start_in) {
                        record_constructor.set_externally_start_start_connected();
                    } else if (h == end_in) {
                        record_constructor.set_externally_start_end_connected();
                    }
                    return true;
                });
                graph->follow_edges(end_out, false, [&](const handle_t& h) {
                    if (h == end_in) {
                        record_constructor.set_externally_end_end_connected();
                    } else if (h == start_in) {
                        record_constructor.set_externally_start_end_connected();
                    }
                    return true;
                });
            }

        }
        first_component_num += temp_index->components.size();
    }

    /* Now go through everything again and give everything children */
    for (size_t temp_index_i = 0 ; temp_index_i < temporary_indexes.size() ; temp_index_i++) {
        TemporaryDistanceIndex* temp_index = temporary_indexes[temp_index_i];
        for (size_t temp_snarl_i = 0 ; temp_snarl_i < temp_index->temp_snarl_records.size() ; temp_snarl_i ++) {
            temp_record_ref_t temp_snarl_ref = make_pair(TEMP_SNARL, temp_snarl_i);
            //Get the temporary index for this snarl
            const TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record = temp_index->get_snarl(temp_snarl_ref);
            if (!temp_snarl_record.is_trivial && !temp_snarl_record.is_simple) {
                //And a constructor for the permanent record, which we've already created
                SnarlRecordWriter snarl_record_constructor (&snarl_tree_records,
                        record_to_offset[make_pair(temp_index_i, temp_snarl_ref)]);
                //Now add the children and tell the record where to find them
                snarl_record_constructor.set_child_record_pointer(snarl_tree_records->size());
                for (temp_record_ref_t child : temp_snarl_record.children) {
                    snarl_record_constructor.add_child(record_to_offset[make_pair(temp_index_i, child)]);

                    //Check if the child is a tip, and if so set start/end_tip connectivity of parent snarl
                    if (child.first == TEMP_NODE) {
                        auto& temp_node_record = temp_index->get_node(child);
                        if (temp_node_record.is_tip) {
                            if (temp_node_record.distance_left_start != std::numeric_limits<size_t>::max() ||
                                 temp_node_record.distance_right_start != std::numeric_limits<size_t>::max()){
                                snarl_record_constructor.set_start_tip_connected();
                            }
                            if (temp_node_record.distance_left_end != std::numeric_limits<size_t>::max() ||
                                 temp_node_record.distance_right_end != std::numeric_limits<size_t>::max()){
                                snarl_record_constructor.set_end_tip_connected();
                            }
                        }
                    } else {
                        auto& temp_chain_record = temp_index->get_chain(child);
                        if (temp_chain_record.is_tip) {
                            if (temp_chain_record.distance_left_start != std::numeric_limits<size_t>::max() ||
                                 temp_chain_record.distance_right_start != std::numeric_limits<size_t>::max()){
                                snarl_record_constructor.set_start_tip_connected();
                            }
                            if (temp_chain_record.distance_left_end != std::numeric_limits<size_t>::max() ||
                                 temp_chain_record.distance_right_end != std::numeric_limits<size_t>::max()){
                                snarl_record_constructor.set_end_tip_connected();
                            }
                        }
                    }
                }
            }
        }
        // This temporary index is fully copied in, so free it
        *temp_index = TemporaryDistanceIndex();
    }


}



}
