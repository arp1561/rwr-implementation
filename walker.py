import sys
import numpy as np
import networkx as nx
from sklearn.preprocessing import normalize

CONV_THRESHOLD = 0.000001

class Walker:                                                               #Walker class

    def __init__(self, original_ppi):            #Contructor function
        self._build_matrices(original_ppi)          

    def run_exp(self, source, restart_prob, og_prob, node_list=[]):         #RWR Function
        self.restart_prob = restart_prob
        self.og_prob = og_prob

        p_0 = self._set_up_p0(source)                                       #sets prob vectors for the seed list
        diff_norm = 1
    
        p_t = np.copy(p_0)

        while (diff_norm > CONV_THRESHOLD):
            p_t_1 = self._calculate_next_p(p_t, p_0)                        #calculates next value 

            diff_norm = np.linalg.norm(np.subtract(p_t_1, p_t), 1)          #calculates current threshold value 

            p_t = p_t_1

        for node, prob in self._generate_rank_list(p_t):
            print '{}\t{:.10f}'.format(node, prob)


    def _generate_rank_list(self, p_t):
        """ Return a rank list, generated from the final probability vector.

        Gene rank list is ordered from highest to lowest probability.
        """
        gene_probs = zip(self.OG.nodes(), p_t.tolist())
        for s in sorted(gene_probs, key=lambda x: x[1], reverse=True):
            yield s[0], s[1]


    def _calculate_next_p(self, p_t, p_0):          
        """ Calculate the next probability vector. """
        epsilon = np.squeeze(np.asarray(np.dot(self.og_matrix, p_t)))
        no_restart = epsilon * (1 - self.restart_prob)
        restart = p_0 * self.restart_prob
        return np.add(no_restart, restart)


    def _set_up_p0(self, source):
        """ Set up and return the 0th probability vector. """
        p_0 = [0] * self.OG.number_of_nodes()
        for source_id in source:
            try:
                source_index = self.OG.nodes().index(source_id)     
                p_0[source_index] = 1 / float(len(source))      #generates the probability
            except ValueError:
                sys.exit("Source node {} is not in original graph. Source: {}. Exiting.".format(
                          source_id, source))
        return np.array(p_0)


    def _build_matrices(self, original_ppi):
        original_graph = self._build_og(original_ppi)


        self.OG = original_graph
        og_not_normalized = nx.to_numpy_matrix(original_graph)
        self.og_matrix = self._normalize_cols(og_not_normalized)
        self.tsg_matrix = None




    def _build_og(self, original_ppi):
        """ Build the original graph, without any nodes removed. """

        try:
            graph_fp = open(original_ppi, 'r')
        except IOError:
            sys.exit("Could not open file: {}".format(original_ppi))

        G = nx.Graph()
        edge_list = []

        # parse network input
        for line in graph_fp.readlines():
            split_line = line.rstrip().split('\t')
            '''
            if len(split_line) > 3:
                # assume input graph is in the form of HIPPIE network
                edge_list.append((split_line[1], split_line[3],
                                  float(split_line[4])))
            '''
            if len(split_line) < 3:
                # assume input graph is a simple edgelist without weights
                edge_list.append((split_line[0], split_line[1], float(1)))
            else:
                # assume input graph is a simple edgelist with weights
                edge_list.append((split_line[0], split_line[1],
                                  float(split_line[2])))

        G.add_weighted_edges_from(edge_list)
        graph_fp.close()

        return G


    def _normalize_cols(self, matrix):
        """ Normalize the columns of the adjacency matrix """
        return normalize(matrix, norm='l1', axis=0)

'''

main function and reading files below

'''
def generate_seed_list(path):
    seed_list = []

    try:
        fp = open(path,"r")
    except IOError:
        sys.exit("Error Opening Seed file! Incorrect path or name")

    for line in fp.readlines():
        info = line.rstrip().split()
        if len(info)>1:
            seed_list.append(info[1])
        else:
            seed_list.append(info[0])
    fp.close()
    print "Successfully opened file and transfered data"
    return seed_list

def main(argv):                                                        #Main Function 

    input_graph = "testdata/test_network.ppi"                          #ENTER Network path here 
    seed_list_path = "testdata/test_seed.txt"                          #Enter SEED LIST PATH HERE     
    seed_list = generate_seed_list(seed_list_path)                     #Opens up seed_file and saves it into seed_list
    
    wk = Walker(input_graph)                      #inits the walker class 
    wk.run_exp(seed_list,restart_prob=0.7,og_prob=0.1,node_list=[])    #Runs RWR




if __name__=='__main__':
    main(sys.argv)
