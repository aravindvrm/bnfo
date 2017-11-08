"""
Aravind Veerappan
BNFO 601 - Exam 2
Question 3. Linkage clustering
"""
from __future__ import division
from math import sqrt


class Tree(object):
    """Defines a full-featured tree class for implementing hierarchical clustering"""

    nodes = {}      # A static attribute for the class.  This dict tracks all of the class instances
                    # keys are in form "GENEX" or "NODEX", while values are instances of Tree (i.e. nodes)..
                    # that is, this class tracks all of its instances by name.
                    # Reminder to look to see if there is a built-in python dict that does this already..

    distance = {}   # Static attribute that tracks the distance between nodes
                    # Keys take the forms of tuples of two node names.
                    # e.g. ("GENE1", GENE20"), ("GENE1", "NODE12"), or ("NODE13", "NODE1")

    dummy = -9999999  # a dummy value for initializing the bestDistance var to something out-of-range of a Pearson

    bestDistance = nextBestDistance = dummy
    bestKeyLeft = bestKeyRight = nextBestKeyLeft = nextBestKeyRight = str
    # Static class attributes tracking the node pair resulting in the minimal distance

    metric = linkage = str
    # Although it is passed as an optional argument to the initializer for user convenience, the instance variable
    # indicating what similarity metric we are using will be re-recorded here so that we may use it in static methods

    # Change linkage type to 'Single' or 'Complete' here
    def __init__(self, internal_name, gene_name, data, metric='Euclidean', right=None,
                 left=None, nodes_represented=None, linkage='Complete'):

        self.internal_name = internal_name
        # This is the name used internally as a key to the nodes dict .. e.g. "GENE2" or "NODE21"
        self.name = gene_name   # This is the gene name from the input file
        self.right = right
        # A key into the nodes dict representing left child -- None if internal name is in form "GENEX"
        self.left = left
        # Same for right child -- None if internal name is in form "GENEX", populated if internal name in form "NODEX"
        self.taken = False   # bool indicating if this node instance has been paired up yet during the clustering phase
                             # New nodes are naturally not yet paired at time of instantiation
        self.data = data     # a list representing the gene expression vector associated with this gene or node
        if not nodes_represented:
            self.nodes_represented = [internal_name]
        else:
            self.nodes_represented = nodes_represented
        # a list with elements corresponding to names/keys of all child nodes of this node instance

        if metric == 'Euclidian':
            self.distance_function = self.__Euclidian_distance
            Tree.metric = 'Euclidian'
        elif metric == 'Pearson':
            self.distance_function = self.__Pearson_correlation
            Tree.metric = 'Pearson'
        elif linkage == 'Single':
            self.distance_function = self.SingleLinkage
            Tree.linkage = 'Single'
        elif linkage == 'Complete':
            self.distance_function = self.CompleteLinkage
            Tree.linkage = 'Complete'
        else:
            self.distance_function = self.__Euclidian_distance
            Tree.metric = 'Euclidian'
            print "Distance metric input was not understood. Defaulting to Euclidian distance"

        # Wow, what did we just do?!? Yes, we can assign functions and methods to variables! In Python, functions,
        # methods, and indeed, classes, are all so-called "first-class objects", and can more-or-less be treated just
        # like any other object!

        print "Instantiating node " + self.internal_name

        # The new instance now calculates the distance between itself and other non-paired instances
        # this could be implemented as a separate private method rather than directly within the constructor

        for other_node in Tree.nodes.values():
            #print "Other node is", other_node
            if other_node.taken is False:  # Check to make sure that we are dealing with an unpaired node

                node_pair = (self.internal_name, other_node.internal_name)
                Tree.distance[node_pair] = self.distance_function(other_node)

                print "Distance between", self.internal_name, "and", other_node.internal_name, "is", \
                    Tree.distance[node_pair]
                Tree.distance[other_node.internal_name, self.internal_name] = Tree.distance[node_pair]
                # The lazy way to deal with reciprocal keys!

                if Tree.distance[node_pair] > Tree.bestDistance:
                # Keep track of the closest pair as we go!!

                        Tree.bestDistance = Tree.distance[node_pair]
                        Tree.bestKeyLeft = self.internal_name
                        Tree.bestKeyRight = other_node.internal_name

        Tree.nodes[self.internal_name] = self   # Finally, add this new instance to the class' static nodes container

    def __Euclidian_distance(self, other_node):          # note that this function is private to the class

            x_vector = self.data
            y_vector = other_node.data

            #print "X_vector", x_vector
            #print "Y_vector", y_vector

            deltas = 0

            for i in xrange(len(x_vector)):  # Assuming that left and right gene expression vectors are equal length

                deltas += (x_vector[i] - y_vector[i])**2

            return -sqrt(deltas)
            # returning the negative of the Euclid only so that remainder of program logic can be left alone.
            # Remember that big Euclidian distance = far apart, while big Pearson = close together

    def __Pearson_correlation(self, other_node):

        x_vector = self.data
        y_vector = other_node.data
        x_sum = 0
        y_sum = 0
        xy_sum = 0
        x_sqr_sum = 0
        y_sqr_sum = 0

        for i in range(len(x_vector)):

            x_sum += x_vector[i]
            y_sum += y_vector[i]
            xy_sum += x_vector[i] * y_vector[i]
            x_sqr_sum += x_vector[i]**2
            y_sqr_sum += y_vector[i]**2

        numerator = len(x_vector) * xy_sum - x_sum * y_sum
        denominator = sqrt(len(x_vector) * x_sqr_sum - x_sum**2) * sqrt(len(x_vector) * y_sqr_sum - y_sum**2)

        r = numerator / denominator

        return r

    @staticmethod
    #  Parameters: the starting node and the .cdt output file
    def DFS(node, cdt_fh):

        data = ""
        if node.left:
            Tree.DFS(Tree.nodes[node.left], cdt_fh)

        if node.right:
            Tree.DFS(Tree.nodes[node.right], cdt_fh)
        # Once a gene is reached, print data to .cdt file
        if node.internal_name.startswith('GENE'):
            for i in node.data[0:]:
                data = str(i) + '\t' + data
            print >> cdt_fh, node.internal_name, '\t', node.internal_name, '\t', node.name, '\t', data

    @staticmethod   # note that this function is private to the class
    def __average(x, y):
        """ Hmmm, this averaging routine just seems to return a simple mean (well, actually, if the numbers
         fed it are log transformed it will correspond to a geometric mean). Either way, this method in no
           way takes into account weights for a node.  Can we fix this?
        """
        return [(x[i] + y[i]) / 2 for i in xrange(len(x))]

    def CompleteLinkage(self, other_node):

        dist = -Tree.dummy
        # iterate through the nodes_represented lists for each node and compare all possible distances between genes
        for gene1 in self.nodes_represented:
            for gene2 in other_node.nodes_represented:
                # find the genes that are farthest from each other, i.e. the lowest correlation
                if Tree.distance[gene1, gene2] < dist:
                    dist = Tree.distance[gene1, gene2]

        return dist

    def SingleLinkage(self, other_node):

        dist = Tree.dummy
        # iterate through the nodes_represented lists for each node and compare all possible distances between genes
        for gene1 in self.nodes_represented:
            for gene2 in other_node.nodes_represented:
                # find the genes that are closest, i.e. the highest correlation
                if Tree.distance[gene1, gene2] > dist:
                    dist = Tree.distance[gene1, gene2]

        return dist

    @staticmethod
    # Parameters: .gtr output file
    def Cluster(gtr_fh):  # This creates a number of composite nodes corresponding to the number of gene nodes less 1

        for i in xrange(len(Tree.nodes) - 1):
            print "\nCreating node", i, \
                "consisting of", Tree.bestKeyLeft, "and", Tree.bestKeyRight, "with distance", Tree.bestDistance

            Tree.nodes[Tree.bestKeyLeft].taken = True
            # First, mark as taken the nodes to the left and right we are about to join
            Tree.nodes[Tree.bestKeyRight].taken = True

            leftData = Tree.nodes[Tree.bestKeyLeft].data
            rightData = Tree.nodes[Tree.bestKeyRight].data
            tempdata = Tree.__average(leftData, rightData)
            print "New vector is", tempdata
            tempname = "NODE" + str(i) + 'X'

            # Keep track of the genes contained in each node
            # Merge the nodes_represented lists when creating a new node
            temp_nodes_rep = Tree.nodes[Tree.bestKeyLeft].nodes_represented \
                + Tree.nodes[Tree.bestKeyRight].nodes_represented
            # print temp_nodes_rep
            if Tree.linkage == 'Single':
                Tree.metric = 'Single'
            else:
                Tree.metric = 'Complete'

            Tree(tempname, tempname, tempdata, Tree.metric, Tree.bestKeyLeft, Tree.bestKeyRight, temp_nodes_rep)
            # Now we must figure out what available pair of genes or nodes are now the closest so that we may join
            # them in the next round!  There are undoubtedly more efficient ways to do this, but this works for now.
            # Print node info to .gtr file
            print >> gtr_fh, tempname, "\t", Tree.nodes[tempname].bestKeyLeft, "\t", Tree.nodes[tempname].bestKeyRight, "\t", Tree.bestDistance

            Tree.bestDistance = Tree.dummy
            for node_pair, cur_distance in Tree.distance.items():

                if cur_distance > Tree.bestDistance:
                    left_one, right_one = node_pair
                    if (Tree.nodes[left_one].taken is False) and (Tree.nodes[right_one].taken is False):
                        # testing to make sure that we are only considering unpaired "available" node pairs
                        Tree.bestDistance = cur_distance
                        Tree.bestKeyLeft = left_one
                        Tree.bestKeyRight = right_one
                        # Reset the class attributes tracking the node pair resulting in the minimal distance

        print "Finished clustering"

        return

    def node_dump(self):
        print "Node name is: " + self.internal_name
        print "Node to left is: " + str(self.left)
        print "Node to right is: " + str(self.right)
        print "Paired node?: " + str(self.taken)
        print "Data associated with node: ", self.data
        print "Nodes this node represents: ", self.nodes_represented
        return


def main():
    """As usual, the main line of the program is pretty compact.  Here I have opted to read in the microrarray
    data file in the main line of the program rather than cluttering up my Tree class.
    """

    i = 0
    # filename = 'ratiodata.txt'
    filename = 'BacillusData2.txt'
    infile = open(filename, 'r')

    # Create the .gtr and .cdt output files
    gtr = open(filename.replace('.txt', '.gtr'), 'w')
    cdt = open(filename.replace('.txt', '.cdt'), 'w')

    # Print the header for the .cdt file
    print >> cdt, 'GID\tUNIQID\tNAME\tWt0A T-1.5\tWt0A T0\tWt0A T2\twtsF T2'

    for line in infile.readlines():
        line = line.strip()
        tempdata = line.split('\t')

        if tempdata[0] != 'GENE':  # Instantiate a new node, but ignore the data label row of the input file

            # Tree("GENE" + str(i), tempdata[0], [float(k) for k in tempdata[1:]], 'Euclidian')
            Tree("GENE" + str(i) + 'X', tempdata[0], [float(k) for k in tempdata[1:]], 'Pearson')
            i += 1

    # for current_node in Tree.nodes:
    #    Tree.nodes[current_node].node_dump()

    Tree.Cluster(gtr)
    Tree.DFS(Tree.nodes['NODE' + str(i-2) + 'X'], cdt)


if __name__ == '__main__':
    main()