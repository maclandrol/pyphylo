/*******************************************************
filename: 	simple_nexus.h
author: 	David Bryant (bryant@math.mcgill.ca)
created: 	August 2002
copyright:	David Bryant 2002 


Routines for readingin the blocks we need from a NEXUS file. This approach is neither
robust nor flexible - but the idea is to have something quick to write and quick to compile.
********************************************************/

#ifndef SIMPLE_NEXUS
#define SIMPLE_NEXUS


#include"global.h"
#include"bit_set.h"

typedef vector<short int> sequence;
enum data_types {standard, DNA, RNA, Protein};
enum transform_types {NONE,HAMMING, GTR, DICE, JAKARD, COVGTR, USER,SITESPEC};


//bool get_word(istream& is, string& next_word, bool make_upper=true);
//bool get_number(istream& is, int& num);
//bool get_double(istream& is, double& num);


short int get_state_id(data_types data_type, char ch);
char get_state(data_types data_type, short int id);
bool state_missing(short int id);

struct bootstrap_params {
    string outfile;
    int replicates;
    bool same_network;
    bool plot_missing;
};

struct st_assumpt_params {
    bool constrained;
    int power;
    double cutoff;
    transform_types transform_type;
    bool equal_rates;
    double shape;
    double pinvar;
    double switch_rate;
    double ss_weight;
    string dist_output;
    string seq_output;
};
    
void write_taxa_block(ostream& os, const vector<string>& taxa_names);
void write_seq_block(ostream& os, data_types seq_type, const vector<string>& taxa_names, const vector<sequence>& seqs);
void simple_bound(data_types seq_type, const vector<string>& taxa_names, const vector<sequence>& seqs);

void write_distances_block(ostream& os, const vector<string>& taxa_names, const dbl_array& D);
void write_st_splits_block(ostream& os, const vector<string>& taxa_names, const vector<bit_set>& splits, const vector<double>& weights);

bool read_nexus_file(	istream& is,
                      string& error_msg,
                      vector<string>& taxa_names,
                      vector<sequence>& seqs,
                      data_types& seq_type,
                      vector<string>& char_names,
                      dbl_array& D,
                      st_assumpt_params& sta_params,
                      bootstrap_params& bs_params);    
#endif



