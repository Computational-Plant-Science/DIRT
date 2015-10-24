'''
kmeans.py

All credits go to to the unknown author who provided the code sniplet.
(http://www.daniweb.com/software-development/python/threads/31449/k-means-clustering)
We just modified tiny pieces.

The code is free for non-commercial use.
Please contact the author for commercial use.

Please cite the DIRT Paper if you use the code for your scientific project.

Bucksch et al., 2014 "Image-based high-throughput field phenotyping of crop roots", Plant Physiology

-------------------------------------------------------------------------------------------
Author: Alexander Bucksch
School of Biology and Interactive computing
Georgia Institute of Technology

Mail: bucksch@gatech.edu
Web: http://www.bucksch.nl
-------------------------------------------------------------------------------------------

Copyright (c) 2014 Alexander Bucksch
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

  * Neither the name of the DIRT Developers nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

'''

'''
# standard python imports
'''
import math, random

# -- The Point class represents points in n-dimensional space
class Point:
    # Instance variables
    # self.coords is a list of coordinates for this Point
    # self.n is the number of dimensions this Point lives in (ie, its space)
    # self.reference is an object bound to this Point
    # Initialize new Points
    def __init__(self, coords, reference=None):
        self.coords = coords
        self.n = len(coords)
        self.reference = reference
    # Return a string representation of this Point
    def __repr__(self):
        return str(self.coords)
    def __getitem__(self,idx):
        return self.coords[idx]
# -- The Cluster class represents clusters of points in n-dimensional space
class Cluster:
    # Instance variables
    # self.points is a list of Points associated with this Cluster
    # self.n is the number of dimensions this Cluster's Points live in
    # self.centroid is the sample mean Point of this Cluster
    def __init__(self, points):
        # We forbid empty Clusters (they don't make mathematical sense!)
        if len(points) == 0: raise Exception("ILLEGAL: EMPTY CLUSTER")
        self.points = points
        self.n = points[0].n
        # We also forbid Clusters containing Points in different spaces
        # Ie, no Clusters with 2D Points and 3D Points
        for p in points:
            if p.n != self.n: raise Exception("ILLEGAL: MULTISPACE CLUSTER")
        # Figure out what the centroid of this Cluster should be
        self.centroid = self.calculateCentroid()
    # Return a string representation of this Cluster
    def __repr__(self):
        return str(self.points)
    def __getitem__(self,idx):
        return self.points[idx]
    # Update function for the <strong class="highlight">K-means</strong> algorithm
    # Assigns a new list of Points to this Cluster, returns centroid difference
    def update(self, points):
        old_centroid = self.centroid
        self.points = points
        self.centroid = self.calculateCentroid()
        return self.getDistance(old_centroid, self.centroid)
    # -- Get the Euclidean distance between two Points
    def getDistance(self,a, b):
        # Forbid measurements between Points in different spaces
        if a.n != b.n: raise Exception("ILLEGAL: NON-COMPARABLE POINTS")
        # Euclidean distance between a and b is sqrt(sum((a[i]-b[i])^2) for all i)
        ret = 0.0
        for i in range(a.n):
            ret = ret + pow((a.coords[i] - b.coords[i]), 2)
        return math.sqrt(ret)
    # Calculates the centroid Point - the centroid is the sample mean Point
    # (in plain English, the average of all the Points in the Cluster)
    def calculateCentroid(self):
        centroid_coords = []
        # For each coordinate:
        for i in range(self.n):
            # Take the average across all Points
            centroid_coords.append(0.0)
            for p in self.points:
                centroid_coords[i] = centroid_coords[i] + p.coords[i]
            if len(self.points)>0: centroid_coords[i] = centroid_coords[i] / len(self.points)
            else: centroid_coords[i] = -1
        # Return a Point object using the average coordinates
        return Point(centroid_coords)
# -- Return Clusters of Points formed by <strong class="highlight">K-means</strong> <strong class="highlight">clustering</strong>
class kMeans:
    def __init__(self,pts):
        self.__points=[]
        for i in pts:
            self.__points.append(Point(i))
            
    def kmeans(self,k, cutoff):
        # Randomly sample k Points from the points list, build Clusters around them
        initial = random.sample(self.__points, k)
        clusters = []
        for p in initial: clusters.append(Cluster([p]))
        # Enter the program loop
        while True:
            # Make a list for each Cluster
            lists = []
            for _ in clusters: lists.append([])
            # For each Point:
            for p in self.__points:
                # Figure out which Cluster's centroid is the nearest
                smallest_distance = self.getDistance(p, clusters[0].centroid)
                index = 0
                for i in range(len(clusters[1:])):
                    distance = self.getDistance(p, clusters[i + 1].centroid)
                    if distance < smallest_distance:
                        smallest_distance = distance
                        index = i + 1
                # Add this Point to that Cluster's corresponding list
                lists[index].append(p)
            # Update each Cluster with the corresponding list
            # Record the biggest centroid shift for any Cluster
            biggest_shift = 0.0
            for i in range(len(clusters)):
                shift = clusters[i].update(lists[i])
                biggest_shift = max(biggest_shift, shift)
            # If the biggest centroid shift is less than the cutoff, stop
            if biggest_shift < cutoff: break
        # Return the list of Clusters
        return clusters
    # -- Get the Euclidean distance between two Points
    def getDistance(self,a, b):
        # Forbid measurements between Points in different spaces
        if a.n != b.n: raise Exception("ILLEGAL: NON-COMPARABLE POINTS")
        # Euclidean distance between a and b is sqrt(sum((a[i]-b[i])^2) for all i)
        ret = 0.0
        for i in range(a.n):
            ret = ret + pow((a.coords[i] - b.coords[i]), 2)
        return math.sqrt(ret)
    # -- Create a random Point in n-dimensional space
    def makeRandomPoint(self,n, lower, upper):
        coords = []
        for _ in range(n): coords.append(random.uniform(lower, upper))
        return Point(coords)

