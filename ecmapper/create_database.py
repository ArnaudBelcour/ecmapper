import urllib 
import os

def download_database(database_folder):
    if not os.path.exists(database_folder):
        os.mkdir(database_folder)

    print('Download Expasy enzyme file')
    urllib.request.urlretrieve('ftp://ftp.expasy.org/databases/enzyme/enzyme.dat', database_folder + '/enzyme.dat')

    print('Download BIGG reactions file')
    urllib.request.urlretrieve('http://bigg.ucsd.edu/api/v2/database_version', database_folder + '/bigg_version.json')
    urllib.request.urlretrieve('http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt', database_folder + '/bigg_models_reactions.txt')

    print('Download ModelSEED reactions file')
    urllib.request.urlretrieve('https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/reactions.tsv', database_folder + '/reactions.tsv')



