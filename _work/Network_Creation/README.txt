Hello, this is the README file for the zebrafish sensory motor pathway (SMP) 
graph theoretic model.  It gives short descriptions of the files needed to make 
a graph representation of the zebrafish SMP.  The functions are written in 
Matlab, which is needed to run the functions, also the software Graphviz, by 
AT&T Bell labs, will need to be installed (nearly any version).

Basically, a giant text file is created that contains the names of neurons and 
their connections to other neurons.  It is a very large text file.  This is then
passed to Graphviz, which lays out the graph and returns an adjacency matrix.  
This adjacency matrix, coupled with a set of names, is the model.

This code was written by Joshua Peterson and Michael Stobb, except where credit
is given to other, smarter people.

###############################################################################

The basics of the system:

Neurons

Number of neurons:
     -Neuromasts: 30 (15 per side)
     -PLL Sensory cells: 40 (20 per side)
     -RS Descending: 290 (145 per side of the RS)
     -Brain interneurons: 906 (3 per Descending cell body, plus a few extra)
     -Spinal Interneurons: 700 (14 per hemisegment for 25 segments)
     -Spinal Motor Neurons: 650 (3 primary and 10 secondary per hemisegment for 25 segments)


Connections

Total PLL connections: 70
All connections are from neuromasts to sensory cells, 35 connections per side.


Total PLL to RS connections: 1,000
Specifications within:
     -Sensory-Brain interneuron: 566
     -Sensory-RS descending: 434


Total RS connections: 70,000
Specifications within the RS:
     -Brain interneuron-Brain interneuron: 35,700
     -RS descending-RS descending: 1,750
     -RS descending-brain interneurons: 32,550


Total RS to Spine connections: 34,800
Specifications within:
     -RS Descending-spinal interneuron: 18,045
     -RS descending-spinal motor neuron: 16,755


Total Spinal connections: 61,244
Specifications within:
     -Spinal interneuron-spinal motor neuron: 43,352
     -Spinal interneuron-Spinal interneuron: 17,892


So the total number of connections is 167,114, and the total number of neurons is 2,616.

###############################################################################

Individual Functions:

activatedCells.mat
	A multi-dimensional list of the neurons shown to be activated during an escape response.

activeCellsList.mat
	A true list of the same information above.

dendriteMap.mat
	A multi-dimensional list of the descending neurons and where their dendrites have been shown to be.  Also included in this list is descending neuron axon length and projection patterns (contralateral, ipsilateral, etc.).

DesNeuronList.mat
	A simple list of all the names of descending neurons, split up by left and right sides.

drawDot.m
	A function written by Bjorn Wielens that loads a Graphviz Dot file and returns an adjacency matrix in Matlab.  Closely examine the help for drawDot to run properly.

lateral_line_stoch.dot
	A simple example setup for the Posterior Lateral Line (PLL) of the zebrafish.  It is later modified.

LLtoRS.m
	A function that writes a Dot file for the connections between the PLL and the Brain of the zebrafish.

makeClusters.m
	A function that makes clusters in Dot files.

makeRSdot2.dot
	A function that creates the zebrafish brain as a Dot file. 

makeSMP.m
	The main function that calls all other necessary functions and outputs an adjacency matrix and list of names corresponding to that matrix.

makeSpinedot.m
	A function that creates the spine of the zebrafish as a Dot file.

RS2Spine3.m
	A function that writes a Dot file for the connections between the Brain and the Spine of the zebrafish.

rsMap.mat
	A multi-dimensional list of the sections of the zebrafish brain and which descending neurons are contained in each.

SpineList.mat
	A simple list of all the types of neurons present in each segment of the spine.

updateMatrixDed.m
	A function that takes the adjacency matrix from drawDot and updates it accordingly.

#############################################################################

To run:

Have all files in the same directory and run the following command in the Matlab command window:

[Adj Names] = makeSMP(file);

Where file is a string of the file name you wish to save the Dot file as, eg file = 'SMPtest.dot', Adj is adjacency matrix output and Names is the list of cell names corresponding to that particular matrix.  This command will probably take several minutes to execute, so you will have to be patient.  Again, it will only work correctly if Graphviz is installed on your system.


