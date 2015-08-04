/*******************************************************
filename: 	simple_nexus.cpp
author: 	David Bryant (bryant@math.mcgill.ca)
created: 	August 2002
copyright:	David Bryant 2002 


Routines for readingin the blocks we need from a NEXUS file. This approach is neither
robust nor flexible - but the idea is to have something quick to write and quick to compile.
********************************************************/

#include"simple_nexus.h"
#include <string.h>


const char* DNA_STATES = "ACGT";
const char* RNA_STATES = "ACUT";
const char* PROTEIN_STATES = "ARNDCEQGHILKMFPSTWYVX";
const char* STANDARD_STATES = "01";
char MISSING = '?';
char GAP = '-';
const short int MISSING_ID = -1;
const short int GAP_ID = -2;
const short int BAD_STATE = -100;


/* Returns true if ch is punctuation 

This code was taken from Paul Lewis's NCL class library*/
static bool is_punctuation(char ch) {
    char punctuation[21];
	punctuation[0]  = '(';
	punctuation[1]  = ')';
	punctuation[2]  = '[';
	punctuation[3]  = ']';
	punctuation[4]  = '{';
	punctuation[5]  = '}';
	punctuation[6]  = '/';
	punctuation[7]  = '\\';
	punctuation[8]  = ',';
	punctuation[9]  = ';';
	punctuation[10] = ':';
	punctuation[11] = '=';
	punctuation[12] = '*';
	punctuation[13] = '\'';
	punctuation[14] = '"';
	punctuation[15] = '`';
	punctuation[16] = '+';
	punctuation[17] = '-';
	punctuation[18] = '<';
	punctuation[19] = '>';
	punctuation[20] = '\0';
    
    return (strchr(punctuation,ch) != NULL);
    
}

static bool is_whitespace(char ch) {
    char whitespace[5];
	whitespace[0]  = ' ';
	whitespace[1]  = '\t';
	whitespace[2]  = '\n';
        whitespace[3]  = '\r';
	whitespace[4]  = '\0';

    return (strchr(whitespace,ch) != NULL);
}

/***********
get_char

Reads in whitespace, then a single character

Function returns false if the EOF is reached before a word is read in.
***************/

static bool get_char(istream& is, char& ch, bool make_upper=true) {
    bool end_of_whitespace;

    /* Flush any whitespace before the word */
    end_of_whitespace = false;

    while (!end_of_whitespace) {
        is.get(ch);
        if ((ch == EOF)||!is) {
            return false;
        }
        else if (is.bad()) {
            return false;
        }
        else if (ch=='[') {
            /* Comment - skip everything until we get to the end of the comment */
            while((is) && (ch!=']'))
                is.get(ch);
            if (!is) {
                return false;
            }
        }
        else if (!is_whitespace(ch))
            end_of_whitespace = true;
    }

    if (make_upper)
        ch=toupper(ch);
    return ch;
}



    
/***********
get_word

Reads in whitespace, then either a single puctuation character
or non-punctuation characters, terminating with the next punctuation character.

Function returns false if the EOF is reached before a word is read in.
***************/

static bool get_word(istream& is, string& next_word, bool make_upper=true) {
    
    char ch;
    bool end_of_whitespace;
    
    /* Flush any whitespace before the word */
    end_of_whitespace = false;
    
    while (!end_of_whitespace) {
        is.get(ch);
	//	fprintf(stderr, "char %c\n",ch);
        if ((ch == EOF)||!is) {
            next_word = "EOF";
            return false;
        }
        else if (is.bad()) {
            next_word = "BADSTREAM";
            return false;
        }
        else if (ch=='[') {
            /* Comment - skip everything until we get to the end of the comment */
            while((is) && (ch!=']'))
                is.get(ch);
            if (!is) {
                next_word = "OPENCOMMENT";
                return false;
            }
        }
        else if (!is_whitespace(ch))
            end_of_whitespace = true;
    }
        
    next_word="";
    if (ch=='\'') { /* Single quote. The word includes everything until the next single quote (but not the quotations!)*/
        while(is) {
            is.get(ch);
            if (ch=='\'')
                return true;
            next_word+=ch;
        }
        next_word="EOF";
        return false;
    }

    
    if (make_upper)
        next_word+=toupper(ch);
    else
        next_word+=ch;
        
   
    if (is_punctuation(ch))  /* Next symbol is punctuation - read in only this character */
        return true;
    
    
        
    while(is) {
        is.get(ch);
	//	fprintf(stderr, "%c\n",ch);
        if (is_punctuation(ch) || is_whitespace(ch)) { /* End of word */
            is.putback(ch);
	    
	    return true;
        }
        else {
            if (make_upper)
                next_word+=toupper(ch);
            else
                next_word+=ch;
        }
    }
    
    next_word = "UNKNOWN";
   
    return false;
}

static bool get_number(istream& is, int& num) {
    string next_word;
    if (!get_word(is,next_word)) return false;
    num = (int) atol(next_word.c_str());
    return true;
}

static bool get_double(istream& is, double& num) {
    string next_word;
    if (!get_word(is,next_word)) return false;
    num =  atof(next_word.c_str());
    return true;
}

static bool skip_to_keyword(istream& is, const string& keyword) {
    
    string next_word="";
    
    while (next_word != keyword)
        if (!get_word(is,next_word))
            return false; /* Problems getting next word */

    return true;
}


/***********
get_next_block name

Skips through everything until a BEGIN keyword is found, then puts the following word into block_name
and parses the following semicolon. If it can't find the block, or runs into other problems, returns false 
(otherwise returns true).

******************/

static bool get_next_block_name(istream& is, bool& no_more_blocks, string& block_name) {
    
    string next_word;
    no_more_blocks = false;
    
    if (!skip_to_keyword(is,"BEGIN")) {
        no_more_blocks = true;
        return true;
    }
    if (!get_word(is,block_name)) return false;
    if (!get_word(is,next_word) || next_word!=";") return false;
    
    return true;
}


/*************
read_taxa_block


Crude parser of a taxa block object.
**************/
static bool read_taxa_block(istream& is, vector<string>& taxa_names) {

    string next_word;
    int i,ntax;
    
    // fprintf(stderr, "read taxa block 1\n");
    if (!get_word(is,next_word) || (next_word != "DIMENSIONS")) return false;
    //fprintf(stderr, "read taxa block 2\n");
    if (!get_word(is,next_word) || (next_word != "NTAX")) return false;
    //fprintf(stderr, "read taxa block 3\n");
    if (!get_word(is,next_word) || (next_word != "=")) return false;
    // fprintf(stderr, "read taxa block 4\n");
    if (!get_number(is,ntax) || (ntax <= 0) ) return false;
    // fprintf(stderr, "read taxa block 5\n");
    if (!get_word(is,next_word) || (next_word != ";")) return false;
    // fprintf(stderr, "read taxa block 6\n");
    if (!get_word(is,next_word) || (next_word != "TAXLABELS")) return false;
    taxa_names.resize(ntax);
    //fprintf(stderr, "read taxa block 7 ntax %i\n",ntax);
    for(i=0;i<ntax;i++) {
     
        if (!get_word(is,taxa_names[i],false))
            return false;
	//fprintf(stderr, "done word %i!\n ",i);
	//fprintf(stderr, "taxa name %s\n",taxa_names[i]);
    }
    // fprintf(stderr, "read taxa block 8\n");
  
    if (!get_word(is,next_word) || (next_word != ";")) return false;
    //   fprintf(stderr, "read taxa block 9\n");
    if (!get_word(is,next_word) || (next_word != "END")) return false;
    // fprintf(stderr, "read taxa block 10\n");
    if (!get_word(is,next_word) || (next_word != ";")) return false;

    return true;    
}

static bool read_dist_block(istream& is, const vector<string>& taxa_names, dbl_array& D) {
    
    string next_word;
    int triangle;
    bool diagonal = true;
    bool labels = true;
    int i,j;
    int ntax;
    
    if (!get_word(is,next_word)) return false;
    if (next_word == "DIMENSIONS") {
        if (!skip_to_keyword(is,";"))
            return false;
        if (!get_word(is,next_word)) return false;
    }
    if (next_word == "FORMAT") {
        for( ; ; ) {
            if (!get_word(is,next_word)) return false;
            if (next_word == ";")
                break;
            else  if (next_word == "TRIANGLE") {
                if (!get_word(is,next_word) || (next_word!="=")) return false;
                if (!get_word(is,next_word)) return false;
                if (next_word == "LOWER")
                    triangle=0;
                else if (next_word == "UPPER")
                    triangle=1;
                else if (next_word == "BOTH")
                    triangle=2;
                else return false;
            }
            else if (next_word=="NO") {
                if (!get_word(is,next_word) ) return false;
                if (next_word == "DIAGONAL")
                    diagonal = false;
                else if (next_word == "LABELS")
                    labels = false;
            }
            else if (next_word=="DIAGONAL") {
                diagonal = true;
            }
	    else if (next_word=="NODIAGONAL") {
	      diagonal = false;
	    }
            else if (next_word=="LABELS") {
                labels = true;
            }
	    else if (next_word=="NOLABELS") {
	      labels = false;
	    }
            else return false;
        }
        
        
        if (!get_word(is,next_word)) return false;
        
    }
    if (next_word!="MATRIX") return false;
    //fprintf(stderr, "WE got to the step with matrix in it.\n");
    ntax = taxa_names.size();
    D.resize(ntax);
    for (i=0;i<ntax;i++) 
        D[i].resize(ntax);
    
    for (i=0;i<ntax;i++) {
        if (labels)
            if (!get_word(is,next_word,false) || next_word!=taxa_names[i]){
                fprintf(stderr, "This throw an error if the name is not found\n" );
                return false;
            }
        if (triangle==0 || triangle == 2) {
            for(j=0;j<i;j++) {
                if(!get_double(is,D[i][j])){
                    fprintf(stderr, "This throw an error if we cannot convert to double\n" );
                    return false;
                } 
                D[j][i] = D[i][j]; //cout<<D[i][j]<<" * ";
            }
        }
        if (triangle==2 || diagonal) {
            if(!get_double(is,D[i][i])){
                fprintf(stderr, "This throw an error if we cannot convert to double wih diagonal enabled\n" );
                return false;
             }
            //cout<<D[i][i]<<" & ";

        }
        else
            D[i][i] = 0.0;
        if (triangle==1 || triangle == 2) {
            for(j=i+1;j<ntax;j++) {
                if(!get_double(is,D[i][j])) return false;
                D[j][i] = D[i][j]; //cout<<D[i][j]<<" ! ";

            }
        }
    }
        
    if (!get_word(is,next_word) || (next_word != ";")) return false;
    if (!get_word(is,next_word) || (next_word != "END")) return false;
    if (!get_word(is,next_word) || (next_word != ";")) return false;
    
    return true;
    
    
}



void write_distances_block(ostream& os, const vector<string>& taxa_names, const dbl_array& D) {

    int i,j;
    int ntax = taxa_names.size();

    os<<"BEGIN DISTANCES;\n";
    os<<"\tFORMAT NO LABELS NO DIAGONAL TRIANGLE=LOWER;";
    os<<"\tMATRIX\n";
    for(i=0;i<ntax;i++) {
        for(j=0;j<i;j++)
            os<<D[i][j]<<" ";
        os<<"\n";
    }
    os<<"\t;\n";

    os<<"END;";

}


 void write_taxa_block(ostream& os, const vector<string>& taxa_names) {
    int i;
    
    int ntax = taxa_names.size();
    os<<"BEGIN TAXA;\n";
    os<<"\tDIMENSIONS NTAX = "<<ntax<<";\n";
    os<<"\tTAXLABELS\n";
    for(i=0;i<ntax;i++) 
        os<<"\t\t "<<taxa_names[i]<<"\n";
    os<<"\t;\n";
    os<<"END;\n\n";
}

void write_st_splits_block(ostream& os, const vector<string>& taxa_names, const vector<bit_set>& splits, const vector<double>& weights) {

    int i,j;
    int ntax = taxa_names.size();
    int nsplits = splits.size();
    
    os<<"BEGIN ST_SPLITS;\n";
    os<<"\tDIMENSIONS NTAX = "<<ntax<<" NSPLITS = "<<nsplits<<";\n";
    os<<"\tFORMAT NO LABELS WEIGHTS;\n";
    os<<"\tPROPERTIES WEAKLY COMPATIBLE CYCLIC;\n";
    os<<"\tMATRIX\n";
    
    for (i=0;i<nsplits;i++) {
        os<<"\t"<<weights[i]<<"\t";
        if (splits[i].member(0)) {
            for (j=0;j<ntax;j++) {
                if (!splits[i].member(j))
                    os<<" "<<j+1;
            }
        }
        else {
            for (j=0;j<ntax;j++) {
                if (splits[i].member(j))
                    os<<" "<<j+1;
            }
        }
        os<<",\n";
    }
    os<<"\t;\n";
    os<<"END;";

}

static bool read_st_assumptions_block(istream& is, st_assumpt_params& params) {
                              
    string next_word;
    

    for( ; ; ) {
        if (!get_word(is,next_word)) return false;
        if (next_word=="END") {
            if (!get_word(is,next_word) || (next_word != ";")) return false;
            return true;
        }
        else if (next_word=="EDGE_TRANSFORM") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_word(is,next_word)) return false;
            if (next_word == "FM_POW1") {
                params.constrained = false;
                params.power = 1;
            }
            else if (next_word == "FM_POW2") {
                params.constrained = false;
                params.power = 2;
            }
            else if (next_word == "LEAST_SQUARES") {
                params.constrained = false;
                params.power = 0;
            }
            else if  (next_word == "NNEG_FM_POW1") {
                params.constrained = true;
                params.power = 1;
            }
            else if (next_word == "NNEG_FM_POW2") {
                params.constrained = true;
                params.power = 2;
            }
            else if (next_word == "NNEG_LEAST_SQUARES") {
                params.constrained = true;
                params.power = 0;
            }
            else return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;

        }
        else if (next_word=="THRESHOLD") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_double(is,params.cutoff)) return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;
        }
        else if (next_word=="DISTANCE") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_word(is,next_word)) return false;
            if (next_word == "HAMMING")
                params.transform_type = HAMMING;
            else if (next_word == "GTR")
                params.transform_type = GTR;
            else if (next_word == "DICE")
                params.transform_type = DICE;
            else if (next_word == "JAKARD")
                params.transform_type = JAKARD;
            else if (next_word == "COVGTR")
                params.transform_type = COVGTR;
            else if (next_word == "SITESPEC")
                params.transform_type = SITESPEC;
            else if ((next_word == "USER") || (next_word == "USERMATRIX") )
                params.transform_type = USER;
            else
                return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;

        }
        else if (next_word=="RATES") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_word(is,next_word)) return false;
            if (next_word == "EQUAL")
                params.equal_rates = true;
            else if (next_word == "GAMMA")
                params.equal_rates = false;
            else
                return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;
        }
        else if (next_word == "SHAPE") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if(!get_double(is,params.shape) || params.shape<0.0) return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;
        }
        else if (next_word == "PINVAR") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if(!get_double(is,params.pinvar) || params.pinvar<0.0) return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;
        }
        else if (next_word == "SWITCHRATE") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if(!get_double(is,params.switch_rate) || params.switch_rate<0.0) return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;
        }
        else if (next_word == "SSWEIGHT") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if(!get_double(is,params.ss_weight) || params.ss_weight<0.0) return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;
        }
        else if (next_word == "DISTOUTPUT") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_word(is,params.dist_output,false)) return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;
        }
        else if (next_word == "SEQOUTPUT") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_word(is,params.seq_output,false)) return false;
            if (!get_word(is,next_word) || (next_word != ";")) return false;
        }        
        else {
            skip_to_keyword(is,";");
        }
    }
}

 short int get_state_id(data_types data_type, char ch) {
    const char *p;
    const char *s;

    if (ch == MISSING)
        return MISSING_ID;
    else if (ch == GAP)
        return GAP_ID;
    
    switch (data_type) {
        case standard: s = STANDARD_STATES; break;
        case DNA: s = DNA_STATES; break;
        case RNA: s = RNA_STATES; break;
        case Protein: s = PROTEIN_STATES; break;
    };

    p = index(s,ch);
    if (p==NULL)
        return BAD_STATE;
    else
        return (long int)p - (long int)s;
}

char get_state(data_types data_type, short int id) {

    if (id == MISSING_ID)
        return MISSING;
    else if (id == GAP_ID)
        return GAP;
    switch (data_type) {
        case standard: return STANDARD_STATES[id]; break;
        case DNA: return DNA_STATES[id]; break;
        case RNA: return RNA_STATES[id]; break;
        case Protein: return PROTEIN_STATES[id]; break;
    };
    return 'x';
}

bool state_missing(short int id) {
    return ((id == MISSING_ID)||(id==GAP_ID));
}



static bool read_char_block(istream& is, const vector<string>& taxa_names, vector<string> & char_names, vector<sequence>& seqs, data_types& data_type) {

    string next_word;
    int i,j,nchar,ntax;
    bool has_labels, labels_left;
    char ch;
    bool transpose;
    
    char_names.clear();
    has_labels = true;
    labels_left = true;
    ntax = taxa_names.size();
    transpose = false;
    data_type = standard;
    
    /* Read the dimensions block */
    if (!get_word(is,next_word)) return false;
    if (next_word != "DIMENSIONS")
        return false;
    if (!get_word(is,next_word)||(next_word != "NCHAR")) return false;
    if (!get_word(is,next_word) || (next_word != "=")) return false;
    
    if (!get_number(is,nchar) || (nchar <= 0) ) return false;
    if (!get_word(is,next_word) || (next_word != ";")) return false;

    /* Read the format block */
    if (!get_word(is,next_word)) return false;
    if (next_word == "FORMAT") {
        for( ; ; ) {
            if (!get_word(is,next_word)) return false;
            if (next_word == ";")
                break;
            else if (next_word == "DATATYPE") {
                if (!get_word(is,next_word) || (next_word != "=")) return false;
                if (!get_word(is,next_word)) return false;
                if (next_word == "STANDARD")
                    data_type = standard;
                else if (next_word == "DNA")
                    data_type = DNA;
                else if (next_word == "RNA")
                    data_type = RNA;
                else if (next_word == "PROTEIN")
                    data_type = Protein;
                else return false;
            }
            else if (next_word=="LABELS") {
                if (!get_word(is,next_word) || (next_word != "=")) return false;
                if (!get_word(is,next_word)) return false;
                if (next_word == "NO")
                    has_labels = false;
                else if (next_word == "LEFT") {
                    has_labels = true;
                    labels_left = true;
                }
                else if (next_word == "RIGHT") {
                    has_labels = true;
                    labels_left = false;
                }
                else
                    return false;
            }
            else if (next_word=="TRANSPOSE")
                transpose = true;
            else if (next_word=="MISSING") {
                if (!get_word(is,next_word) || (next_word != "=")) return false;
                if (!get_char(is,MISSING)) return false;
                if (MISSING == '\'') {
                    if (!get_char(is,MISSING)) return false;
                    if (!get_char(is,ch) || (ch != '\'')) return false;
                }
            }
            else if (next_word=="GAP") {
                if (!get_word(is,next_word) || (next_word != "=")) return false;
                if (!get_char(is,GAP)) return false;
                if (GAP == '\'') {
                    if (!get_char(is,GAP)) return false;
                    if (!get_char(is,ch) || (ch != '\'')) return false;
                }
            }
            else
                return false;
        }
        if (!get_word(is,next_word)) return false;
    }

    if (next_word != "MATRIX")
        return false;

    seqs.resize(ntax); /* Allocate memory for the sequences */
    for(i=0;i<ntax;i++) 
        seqs[i].resize(nchar);
    
    if (!transpose) { /* Read in the matrix untransposed */
        for(i=0;i<ntax;i++) {
            if (has_labels && labels_left)
                if (!get_word(is,next_word,false) || next_word != taxa_names[i] ) {cout<<taxa_names[i]<<endl;return false;}
            for(j=0;j<nchar;j++) {
                if (!get_char(is,ch)) return false;
                seqs[i][j] = get_state_id(data_type,ch);
                if (seqs[i][j]<-2) return false;
            }
            if (has_labels && !labels_left)
                if (!get_word(is,next_word,false) || next_word != taxa_names[i] ) return false;
        }
    }
    else { /* Read in the matrix transposed */
        if (has_labels)
            char_names.resize(nchar);
        for(j=0;j<nchar;j++) {
            if (has_labels && labels_left)
                if (!get_word(is,char_names[j],false) ) return false;
            for(i=0;i<ntax;i++) {
                if (!get_char(is,ch)) return false;
                seqs[i][j] = get_state_id(data_type,ch);
                if (seqs[i][j]<-2) return false;
            }
            if (has_labels && !labels_left)
                if (!get_word(is,char_names[j],false) ) return false;

        }
    }
    
    if (!get_word(is,next_word) || next_word != ";") return false;

    if (!get_word(is,next_word) || next_word != "END") return false;

    if (!get_word(is,next_word) || next_word != ";") return false;

    return true;
}


static bool read_bootstrap_block(istream& is, bootstrap_params& bs_params) {

    string next_word;

    for ( ; ; ) {
        if (!get_word(is,next_word)) return false;
        if (next_word == "OUTFILE") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_word(is,bs_params.outfile,false)) return false;
        }
        else if (next_word == "REPLICATES") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if ((!get_number(is,bs_params.replicates))||(bs_params.replicates<1)) return false;
        }
        else if (next_word == "SAME") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_word(is,next_word)) return false;
            if (next_word == "YES")
                bs_params.same_network = true;
            else if (next_word == "NO")
                bs_params.same_network = false;
            else
                return false;
        }
        else if (next_word == "PLOTMISSING") {
            if (!get_word(is,next_word) || (next_word != "=")) return false;
            if (!get_word(is,next_word)) return false;
            if (next_word == "YES")
                bs_params.same_network = true;
            else if (next_word == "NO")
                bs_params.same_network = false;
            else
                return false;
        }
        else if (next_word == "END")
            break;
        else
            return false;

        if (!get_word(is,next_word) || (next_word != ";")) return false;
    }
    if (!get_word(is,next_word) || (next_word != ";")) return false;

    return true;
}
        

bool read_nexus_file(	istream& is,
                      string& error_msg,
                      vector<string>& taxa_names,
                      vector<sequence>& seqs,
                      data_types& seq_type,
                      vector<string>& char_names,
                      dbl_array& D,
                      st_assumpt_params& sta_params,
                      bootstrap_params& bs_params) {
    
    string next_word;
    bool found_taxa = false;
    bool found_distances = false;
    bool found_characters = false;
    bool no_more_blocks;
    seqs.clear();
    D.clear();
    
    if (!get_word(is,next_word)) {
        error_msg = "Could not read from the file";
        return false;
    }
    if (next_word!="#NEXUS") {
        error_msg = "The file does not begin with #NEXUS";
        return false;
    }

    for ( ; ; ) {
        if (!get_next_block_name(is,no_more_blocks,next_word)) {
            error_msg = "Error reading header of nexus block";
            return false;
        } 

        if (no_more_blocks) {
	  if (!found_distances && !found_characters) {
	    error_msg = "Could not find or read a DISTANCES or CHARACTERS block in the input file";
	    return false;
	  }
	  return true;
        }
  
        if (next_word=="TAXA") {
            if (!read_taxa_block(is,taxa_names)) {
                error_msg = "Error reading taxa block";
		fprintf(stderr, "Error reading TAXA block.  Please ensure that the taxa names contain no punctutation (for example ':' or '.')\n");
                return false;
            }
            found_taxa = true;
        }
        else if (next_word == "DISTANCES") {
            if (!found_taxa) {
                error_msg = "Error: the TAXA block must appear before the DISTANCES block";
                return false;
            }
            if (!read_dist_block(is,taxa_names,D)) {
                error_msg = "Error reading DISTANCES block";
                return false;
            }
	    found_distances = true;
        }
        else if (next_word == "ST_ASSUMPTIONS") {
            if (!read_st_assumptions_block(is,sta_params)) {
                error_msg = "Error reading ST_ASSUMPTIONS block";
                return false;
            }
        }
        else if (next_word == "CHARACTERS") {
            if (!found_taxa) {
                error_msg = "Error: the TAXA block must appear before the DISTANCES block";
                return false;
            }
            if (!read_char_block(is,taxa_names,char_names, seqs,seq_type)) {
                error_msg = "Error reading CHARACTERS block";
                return false;
            }
            found_characters = true;
        }
        else if (next_word == "BOOTSTRAP") {
            if (!read_bootstrap_block(is,bs_params)) {
                error_msg = "Error reading BOOTSTRAP block";
                return false;
            }
        }
        else {
            if (!skip_to_keyword(is,"END") || !get_word(is,next_word) || (next_word != ";")) {
                error_msg = "Error reading ";
                error_msg += next_word;
                error_msg += " block";
                return false;
            }
        }    
    }
    

}
        


