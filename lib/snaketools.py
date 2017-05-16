#!/usr/bin/env python3
# -*- coding:utf8 -*-

"""
input: directory containing fastq files
output: 
    dictionnary of sorted files 
    file samples.yml (with -c parameters)

samples.yml looks like : 
---------------------------------------
fastq_dir:
    input/raw-data
samples:
    ssr1:
        assay: total
        specie: human
        paired: True
        stranded: True
    ssr2:
        assay: polyA
        specie: human
        paired: True
        stranded:   False
exluded:
    file.fastq:
        reason: 'bad extension (not fastq.gz)ile'
    file.other
        reason: 'bad extension (not fastq.gz)'
    file.fastq.gz
        reason: 'file corrupted'
---------------------------------------
"""

__appname__ = "snaketools"
__shortdesc__ = "bunch of tools for snakemake"
__licence__ = "GPL"
__version__ = "0.4"
__author__ = "Benoit Guibert <benoit.guibert@free.fr>"

import os
import argparse
import yaml

class BuildSamplesFile:
    """ 
    Find fastq files in a directory, and build sample.yml for snakemake and 
    return sample list 
    """
    
    def __init__ (self, fastq_dir, output_file = 'samples.yml'):
        """ Class initialiser """
        self.fastq_dir = fastq_dir                                  # fastq directory
        self.output_file = output_file                              # output file        
        self.samples_dict = {}                                      # initialyze sample dictionary
        self.__add_fastq_dir__()                                    # append fastq_dir
        self.__add_samples__(self.samples_dict)                     # append samples
        self.__set_species__()                                      # None, 'human', 'mouse', ...
        self.__set_stranded__()                                     # None, True, False
        self.__set_assay__()                                        # None, 'polyA', 'total'
        self.__set_file_yaml__(self.output_file, self.samples_dict) # create yaml sample file for snakemake
        

    def __add_fastq_dir__ (self):
        """ Append fastq directory in sample dict """
        self.samples_dict = { 'fastq_dir': self.fastq_dir }


    def __add_samples__ (self, samples_dict):
        """ Append samples to samples dictionnary """
        samples_dict['samples'] = {}
        # eliminate files with bad extension
        candidate_files = self.__exclude_non_fastq_ext__(samples_dict)
        # sort paired samples
        self.__sort_paired_samples__(candidate_files, samples_dict)


    def __exclude_non_fastq_ext__(self, samples_dict):
        # list all files in directory
        file_list = os.listdir(samples_dict['fastq_dir'])
        candidate_files = []
        for file in file_list: 
            split_file = file.split('.')
            if (len(split_file) < 3) or (split_file[-1] != 'gz') or (split_file[-2] != 'fastq') :
                # exclude files with bad extension
                if not 'excluded' in samples_dict:
                    samples_dict['excluded'] = {}
                samples_dict['excluded'][file] = { 'reason': 'bad extension (not fastq.gz)' }
            else:
                # append file as sample candidate
                candidate_files.append(file)
        return candidate_files


    def __sort_paired_samples__(self, candidate_files, samples_dict):
        # print(candidate_files)
        for file in candidate_files:
            basename = ".".join(file.split('.')[:-2])
            if not '_' in basename:
                samples_dict['samples'][basename] = { 'paired': False}
            else:
                basename_split = basename.split('_')
                paired_candidate = "_".join(basename_split[:-1])
                if ( basename_split[-1] == '1' ) and ( paired_candidate+'_2.fastq.gz' in candidate_files):
                        if not paired_candidate in samples_dict['samples']:
                            samples_dict['samples'][paired_candidate] = {'paired': True}
                elif ( basename_split[-1] == '2' ) and ( paired_candidate+'_1.fastq.gz' in candidate_files):
                        if not paired_candidate in samples_dict['samples']:
                            samples_dict['samples'][paired_candidate] = {'paired': True}
                else:
                    samples_dict['samples'][basename] = {'paired': False}


    def __set_species__ (self):
        """ Values: None, 'Human', 'mouse', etc. """
        for sample in self.samples_dict['samples']:
            self.samples_dict['samples'][sample]['specie'] = None

        
    def __set_stranded__ (self):
        """ Values:  None, True, False"""
        for sample in self.samples_dict['samples']:
            self.samples_dict['samples'][sample]['stranded'] = None
            
    
    def __set_assay__(self):
        """ Values: None, 'polyA', 'total' """
        for sample in self.samples_dict['samples']:
            self.samples_dict['samples'][sample]['assay'] = None
        

    def __set_file_yaml__ (self, output_file, samples):
        """ Write yaml output file"""
        with open(output_file, 'w') as outfile:
            yaml.dump(samples, outfile, default_flow_style=False, allow_unicode=True)


class GetSamplesFromFile:
    """ Class doc """
    
    def __init__ (self, sample_file='samples.yml'):
        """ Class initialiser """
        self.sample_file = sample_file
        # Read smaple file
        self.__getSampleFile__()
    
    
    def __getSampleFile__(self):
        """ Function doc """
        with open(self.sample_file) as stream:
            self.samples_dict = yaml.load(stream)
    

    def get_directory(self):
        return self.samples_dict['fastq_dir']


    def get_all(self):
        return self.samples_dict


    def get_sample(self, sample):
        return self.samples_dict['samples'][sample]


    def get_samples(self):
        return self.samples_dict['samples']


    def get_paired(self):
        samples = []
        for sample in self.samples_dict['samples']:
            if self.samples_dict['samples'][sample]['paired']:
                samples.append(sample)
        return samples


    def get_unpaired(self):
        samples = []
        for sample in self.samples_dict['samples']:
            if not self.samples_dict['samples'][sample]['paired']:
                samples.append(sample)
        return samples
   
    
    def get_stranded(self):
        samples = []
        for sample in self.samples_dict['samples']:
            if self.samples_dict['samples'][sample]['stranded']:
                samples.append(sample)
        return samples
   
    
    def get_unstranded(self):
        samples = []
        for sample in self.samples_dict['samples']:
            if not self.samples_dict['samples'][sample]['stranded']:
                samples.append(sample)
        return samples


    def get_assay(self, sample):
        return self.samples_dict['samples'][sample]['assay']


    def get_specie(self):
        pass

def arguments():
    """ parse arguments, in standalone mode"""
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_dir",
                        help = "fastq directory",
                        )
    parser.add_argument("-o", "--output_file",
                        help = "yaml formatted output file, default : sample.yml",
                        default = 'samples.yml',
                        metavar = ('output_file'),
                        )
    parser.add_argument("--version", "-v",
                        action = "version",
                        version = __appname__ + ": " + __version__ ,
                        help = "show version",
                        )
    return parser.parse_args()


if __name__ == "__main__":
    args = arguments()
    samples = BuildSamplesFile(args.fastq_dir, args.output_file)
    

