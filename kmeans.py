import math, random

# All credits go to to the unknown author who provided the code sniplet.
# (http://www.daniweb.com/software-development/python/threads/31449/k-means-clustering)
# We just modified tiny pieces.


class Point:
    """
    A point in n-dimensional space.
    """

    def __init__(self, coords, reference=None):
        self.coords = coords        # list of coordinates for this Point
        self.n = len(coords)        # number of dimensions this Point lives in (ie, its space)
        self.reference = reference  # object bound to this Point

    def __repr__(self):
        """
        Returns a string representation of this point.
        """
        return str(self.coords)

    def __getitem__(self, idx):
        return self.coords[idx]


def get_euclidean_distance(a: Point, b: Point) -> float:
    """
    Get the Euclidean distance between two Points
    """

    if a.n != b.n:
        raise ValueError("Cannot compare points with different dimensionality")

    # Euclidean distance between a and b is sqrt(sum((a[i]-b[i])^2) for all i)
    ret = 0.0
    for i in range(a.n):
        ret = ret + pow((a.coords[i] - b.coords[i]), 2)

    return math.sqrt(ret)


class Cluster:
    """
    A cluster of points in n-dimensional space.
    """

    def __init__(self, points):
        # forbid empty clusters (they don't make mathematical sense!)
        if len(points) == 0:
            raise ValueError("Clusters must not be empty")

        self.points = points    # points associated with this cluster
        self.n = points[0].n    # number of dimensions this cluster's points live in

        # forbid clusters containing points in different-dimensional spaces
        # (e.g., no clusters with both 2D points and 3D points)
        for p in points:
            if p.n != self.n:
                raise ValueError("All points in a cluster must have the same number of dimensions")

        self.centroid = self.calculate_centroid()    # sample mean point of this cluster

    # Return a string representation of this Cluster
    def __repr__(self):
        return str(self.points)

    def __getitem__(self, idx):
        return self.points[idx]

    # Update function for the <strong class="highlight">K-means</strong> algorithm
    # Assigns a new list of Points to this Cluster, returns centroid difference
    def update(self, points):
        old_centroid = self.centroid
        self.points = points
        self.centroid = self.calculate_centroid()
        return get_euclidean_distance(old_centroid, self.centroid)

    # Calculates the centroid Point - the centroid is the sample mean Point
    # (in plain English, the average of all the Points in the Cluster)
    def calculate_centroid(self):
        centroid_coords = []
        # For each coordinate:
        for i in range(self.n):
            # Take the average across all Points
            centroid_coords.append(0.0)
            for p in self.points:
                centroid_coords[i] = centroid_coords[i] + p.coords[i]
            if len(self.points) > 0:
                centroid_coords[i] = centroid_coords[i] / len(self.points)
            else:
                centroid_coords[i] = -1

        # Return a Point object using the average coordinates
        return Point(centroid_coords)


class KMeans:
    def __init__(self, pts):
        self.__points = []
        for i in pts:
            self.__points.append(Point(i))

    def kmeans(self, k, cutoff):
        print('Do the kmeans :-)')

        # Randomly sample k Points from the points list, build Clusters around them
        initial = random.sample(self.__points, k)
        clusters = []
        for point in initial:
            clusters.append(Cluster([point]))

        while True:
            # Make a list for each cluster
            lists = []
            for _ in clusters:
                lists.append([])

            for point in self.__points:
                # Figure out which Cluster's centroid is the nearest
                smallest_distance = get_euclidean_distance(point, clusters[0].centroid)
                index = 0
                for i in range(len(clusters[1:])):
                    distance = get_euclidean_distance(point, clusters[i + 1].centroid)
                    if distance < smallest_distance:
                        smallest_distance = distance
                        index = i + 1
                # Add this Point to that Cluster's corresponding list
                lists[index].append(point)

            # Update each Cluster with the corresponding list
            # Record the biggest centroid shift for any Cluster
            biggest_shift = 0.0
            for i in range(len(clusters)):
                shift = clusters[i].update(lists[i])
                biggest_shift = max(biggest_shift, shift)

            # If the biggest centroid shift is less than the cutoff, stop
            if biggest_shift < cutoff:
                break

        return clusters

    @staticmethod
    def __make_random_point(n, lower, upper):
        """
        Create a random Point in n-dimensional space
        """

        coords = []
        for _ in range(n): coords.append(random.uniform(lower, upper))
        return Point(coords)
