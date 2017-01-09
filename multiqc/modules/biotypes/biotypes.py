#!/usr/bin/env python

""" MultiQC module to parse table with biotype information """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Biotypes', anchor='biotypes',
        href="",
        info="custom table with biotype counts in one column.")
        
        # Find and load any biotype reports
        self.biotype_data = dict()
        for f in self.find_log_files(config.sp['biotypes']):
            self.parse_biotypes(f)
            self.add_data_source(f)

        if len(self.biotype_data) == 0:
            log.debug("Could not find any biotype data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.biotype_data)))
        self.write_data_file(self.biotype_data, 'multiqc_biotype')        

        # Make barplot
        self.intro += self.biotype_barplot()


        
    def parse_biotypes (self, f):
        """ Parses  a file with 2 columns, first one with names, second one with numbers. """
        s_name = f['s_name']
        s_name = self.clean_s_name(s_name, f['root'])
        if s_name in self.biotype_data:
            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        self.biotype_data[s_name] = dict()    
        for l in f['f'].split("\n"):
            ll = re.split('\t| |;|,',l.strip())
            if len(ll) > 1 and ll[1].isdigit():
                self.biotype_data[s_name][ll[0]]=int(ll[1])
        s_name = None

    def biotype_barplot(self):
        # Make a normalised percentage version of the data
        self.biotype_percent = dict()
        headers = []
        for s_name in self.biotype_data:
            self.biotype_percent[s_name]=dict()
            total = sum( self.biotype_data[s_name].values() )
            for k, v in self.biotype_data[s_name].items():
                self.biotype_percent[s_name][k] = (v/total)*100
                if k not in headers: headers.append(k)

        # if biotypes contains spike-in - make another percentage plot without spike-in
        spike_regexp = ['spike-in','ERCC','ercc','spikein']
        spike_name = list(set(spike_regexp) & set(headers))
        if len(spike_name) ==1:
            non_spike = [k for k in headers if k not in spike_name]
            self.biotype_nonspike_percent=dict()
            for s_name in self.biotype_data:
                self.biotype_nonspike_percent[s_name]=dict()
                total = sum( v for k,v in self.biotype_data[s_name].items() if k in non_spike )
                for k, v in self.biotype_data[s_name].items():
                    if not k in non_spike: continue
                    self.biotype_nonspike_percent[s_name][k] = (v/total)*100
            # also, reorder items so that spike-in comes first
            headers = spike_name + non_spike
            
        # make keys from all headers
        keys = OrderedDict()
        for h in headers:
            keys[h]={'name':h}
        #list datasets to use and make labels
        datasets = [self.biotype_data, self.biotype_percent]
        data_labels = [
            {'name': 'Counts', 'ylab': 'Counts'},
            {'name': 'Percentages', 'ylab': 'Percentage'}
        ]
        use_keys = [keys,keys]
        if len(spike_name)==1:
            datasets.append(self.biotype_nonspike_percent)
            data_labels.append({'name': 'Percentages wo spike-in', 'ylab': 'Percentage'})
            use_keys.append(keys)

        # Config for the plot
        pconfig = {
            'id': 'biotype_distribution_plot',
            'title': 'Biotypes',
            'ylab': 'Counts',
            'xlab': 'Biotype',
            'data_labels': data_labels,
            'cpswitch': False
        }
        return plots.bargraph.plot(datasets, use_keys, pconfig)

        
            
            
                    
                    

