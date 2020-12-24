import csv
from collections import OrderedDict


def any_true(ordered_dict: OrderedDict):
    return any(v for v in ordered_dict.values())


def get_traits(path='traits.csv'):
    with open(path, 'U') as file:
        traits = []
        dialect = csv.Sniffer().sniff(file.read(1024))
        file.seek(0)
        reader = csv.reader(file, dialect=dialect)
        for row in reader:
            if row[0] == 'TRAIT':
                continue
            if bool(int(row[1])):
                traits.append(row[0])

    print(f"Including traits: {traits}")

    return traits


def print_header():
    print('------------------------------------------------------------')
    print('DIRT 1.1 - An automatic highthroughput root phenotyping platform')
    print('(c) 2014 Alexander Bucksch - bucksch@uga.edu')
    print('Web application by Abhiram Das - abhiram.das@gmail.com')
    print(' ')
    print('http://dirt.iplantcollaborative.org')
    print(' ')
    print('University of Georgia')
    print('------------------------------------------------------------')
    print(' ')
    print('Initializing folder structure')