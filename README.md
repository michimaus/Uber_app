# Uber_app
***Using the skiplist in this one***

                  __--__ Duican Mihnea-Ionut 314CA __--__


    
    
First of all, I had to read the input data and to add in the hash_table
the names for the streets, in order to store them in a proper way(making 
eficient queries when looking for a specific street).

After this I built the graph using the std::vector and the skip_list(thought 
would be helpful for the future tasks).
NOTE: I modified the skiplist from "Tema 1 SD" quite a bit, as now I am using 
level pointers in stead of full structured nodes to be placed on a level.


TASK 1 & TASK 2
-- Used breadth first search to find the way from the starting point to the finish

TASK 3
-- The main reason I chose to use the skiplist is to have the posibility to make
eficient updates when talking about the links between the nodes
-- Finishing the task I mapped all the distances in order to be able to make
quairies for the fourth task

TASK 4
-- I used three skipl_lists in order to have sorted at any given time my drivers
accourding to the datasets given in the task. For finding the distances necessary
when picking a ride I checked the data I mapped at the end of the task 3.

TASK 5
-- Once again I used data mapped at the end of thask 3. Given the last positian
the driver is to be found I checked all the possible nodes that can be reached
by the order of the amount of fuel the care is going to have when the driver gets
there. Because of using the skiplist I was keeping all the time the nodes in the
needed order, resulting in the end for the right output.
