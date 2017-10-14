"""
ROSALIND - Rabbits and Recurrence Relations

Given: Positive integers n less than 40 and k less than 5.

Return: The total number of rabbit pairs that will be present after n months if we begin
with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter
of k rabbit pairs (instead of only 1 pair).
"""
from __future__ import print_function
from __future__ import division


# modified fibonacci sequence: Fn = Fn-1 + Fn-2 * k
def rabbits(generations, offspring):
    if generations == 0:
        return 0
    elif generations == 1:
        return 1
    else:
        return (rabbits(generations-1, offspring)
                + rabbits(generations-2, offspring) * offspring)


print(rabbits(31, 3))
