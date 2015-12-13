#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
using namespace std;


//a struct called DNA that holds a vector of strings, a string of DNA and the DNA's GC
struct DNA {
    vector<string> headers;
    string aDNA;
    double GC;  
} ;

//declares the constant size of the array
const int MAX_DNA_SIZE = 1000;


//all of the prototypes to the functions
bool isValidDNAchar(char ACGT);
bool isValidDNAstring(string DNA);
string transcribe(string DNA);
double gc(string DNA);
char complementChar(char base);
string complementStr(string DNA);
vector<string> codons(string RNA);
int numCodon(vector<string> codonV,string codon);


int main(){
    ifstream DNAfile;
    ofstream GCcontent,expandedDNA,headerLinesFreq,badDNA;
    //open the input file
    DNAfile.open("dna.txt");
    //creates the new files in the same directory
    GCcontent.open("GCcontent.txt");
    expandedDNA.open("expandedDNA.txt");
    headerLinesFreq.open("headerLinesFreq.txt");
    badDNA.open("badDNA.txt");

    DNA DNAtemp;
    DNA DNAarray[MAX_DNA_SIZE];
    vector<string> hold3;
    string aString,hold2;
    int i =0,actualSize=1,counter=0,largest,freq,smallestIndex,smallestValue,check;
    double hold1;
    
    //reads every string in the file until it is at the end of the file
    while(!DNAfile.eof()){
        // gets one string in the file and stores it in aString which goes through various conditionals
        getline(DNAfile,aString);
        if(aString[0] == '>'){
            //adds the headers into the headers vector in DNAtemp;
            DNAtemp.headers.push_back(aString);
        }
       
        else if (isValidDNAstring(aString) == true ){
            //if it is a valid DNA sequence, then record the header vector,gc and DNA strand of that sequence into the array, then empty each member of DNAtemp to be used again in finding the next DNA sequence
            DNAtemp.GC = gc(aString);
            DNAtemp.aDNA =aString; 
            DNAarray[i] = DNAtemp;     
            i++;
            actualSize++;
            DNAtemp.GC =0;
            DNAtemp.aDNA ="";
            while (DNAtemp.headers.size() > 0){
                DNAtemp.headers.pop_back();
            }
            
        }

        
        else{
            //if the DNA sequence was invalid, then take that DNA sequence and add it to the badDNA file along with its header names
            while (DNAtemp.headers.size() > 0){
                badDNA << DNAtemp.headers[0];
                DNAtemp.headers.pop_back();
            }
            badDNA << aString;
        }
 
        
        
    }
    i=0;
    
    //finds the maximum header frequency 
    while(i < actualSize){
        if (DNAarray[i].headers.size() > largest){
            largest = DNAarray[i].headers.size();
        }
        i++;
    }
    i=0;
    
    // sorts and prints the number of header lines along with their frequencies 
    for (int j = 1; j <= largest; j++){
        while(i < actualSize){
            //checks how many headers an element has
            if (DNAarray[i].headers.size() == j){
                freq++;
            }
            //moves on to the next element
            i++;
        }
        i=0;
        // prints the element header number along with it's frequency 
        headerLinesFreq << j << "  " << freq << endl;
        freq = 0;
        
    }

    

    // sorts the list 
    
    // starts at the first element then reads it throughout the whole array
    for (int i = 0; i < actualSize; i++){
        
        //assumes the smallestindex is the current index and it contains the smallest values 
        smallestIndex = i;
        hold1 = DNAarray[i].GC;
        hold2 = DNAarray[i].aDNA;
        hold3 = DNAarray[i].headers;
        
        // starts a for loop, starting on the index after the current index, so i+1,up to the last index
        for(int j = i + 1; j < actualSize; j++) {
            // if the current value of DNAarray[j].GC is less than the smallest value
            if(DNAarray[j].GC < DNAarray[smallestIndex].GC ){
                //record that index as the index with the smallest value
                smallestIndex = j;  
            }
            // the for loop above shall check the where the smallest GC is and what index it is at compared to the current index
            // 
        }
        
        // if the value of the smallestIndex is changed
        if (smallestIndex != i ){
            // performs the switches of values
            DNAarray[i].GC = DNAarray[smallestIndex].GC;
            DNAarray[i].aDNA = DNAarray[smallestIndex].aDNA;
            DNAarray[i].headers = DNAarray[smallestIndex].headers;
            
            // the initial values of DNAarray[i] are now in the elements DNAarray[smallestIndex]
            DNAarray[smallestIndex].GC = hold1;
            DNAarray[smallestIndex].aDNA = hold2;
            DNAarray[smallestIndex].headers = hold3;
        }
    }
    // stuffs the values in format in the GCcontent file
    for(int i = 1; i < actualSize;i++){
        for (int j = 0; j < DNAarray[i].headers.size();j++){
            GCcontent << DNAarray[i].headers[j] << endl;
        }
        GCcontent << "> GC content: " << DNAarray[i].GC  << endl;
        GCcontent << DNAarray[i].aDNA << endl;
        GCcontent <<  endl;
    }
    
    // stuffs the values in format in the expandedDNA file
    string RNAtemp;
    for(int i = 1; i < actualSize;i++){
        // prints the headers
        for (int j = 0; j < DNAarray[i].headers.size();j++){
            expandedDNA << DNAarray[i].headers[j] << endl;
        }
        // prints the GC content and the DNA string
        expandedDNA << "> GC content: " << DNAarray[i].GC  << endl;
        expandedDNA << DNAarray[i].aDNA << endl;
        
        // prints the transcribed DNA sequence, the number of UGG codons in the RNA sequence and the complement DNA of the current DNA sequence
        RNAtemp = transcribe(DNAarray[i].aDNA);
        expandedDNA << "> RNA: " << RNAtemp << endl;
        expandedDNA << "> Number of UGG(tryptophan): " << numCodon(codons(RNAtemp),"ugg") << endl;
        expandedDNA << "> DNA complement: " << complementStr(DNAarray[i].aDNA) << endl;
        expandedDNA <<  endl;
    }

    
    
    
    
    
     //closes all the streams
    GCcontent.close();
    expandedDNA.close();
    headerLinesFreq.close();
    badDNA.close();
    DNAfile.close();



    return 0;
}

//if a char argument is one of 'a','c','g' or 't', otherwise the function returns false
bool isValidDNAchar(char ACGT){
    if (ACGT == 'a' || ACGT == 'c' || ACGT == 'g' || ACGT == 't'){
        return true;   
    }
    else {
        return false;   
    }
    
}

//returns true  if a string is only made up of valid DNA characters. If any other character exists in the string false is returned
bool isValidDNAstring(string DNA){
    int length = DNA.length(),counter = 0 ;
    
    //iterates through each char in the string, if it hits an invalid char, or it reaches the end of the string, leave the loop
    while (counter < length && isValidDNAchar(DNA[counter]) == true){
        
    counter ++;   
    }
    
    // if the while loop was succesful and didnt find an invalid character,return true, else return false
    if(counter == length){
        
        return true;   
    }
    else{
 
        return false;   
    }
}
//Takes a DNA string argument, and returns its RNA equivalent.
string transcribe(string DNA){
    int length = DNA.length();
    string RNA = "";
    
    // goes through each char in a string, if a char is a 't',then you place a 'u' in the RNA,else, place the char in RNA. 
    for (int i=0;i < length; i++){
        if(DNA[i] == 't'){
            RNA += 'u';   
        }
        else{
            RNA += DNA[i];   
        }
    }
    return RNA;
}

//return a double between 0 and 1 that is the GC content of a DNA sequence
double gc(string DNA){
    int length = DNA.length();
    double total;
    // goes through the whole string, if the current char is a c or g, add one to total
    for (int i = 0; i < length;i++){
        if (DNA[i] == 'c' || DNA[i] == 'g'){   
            total ++;   
        }
    }
    
    //return the GC content
    return total/length;
    
    
    
}


// returns the complement of the base character
// Precondition: the base character argument is a valid DNA character.
char complementChar(char base){
    char complement;
    // makes sure the base is a valid letter
    // if it is a valid base,place the complement value into complement

    if (base == 'a'){
        complement = 't';   
    }
    else if (base == 'c'){
        complement = 'g';   
    }
    else if (base == 'g'){
        complement = 'c';   
    }
    else if (base == 't'){
        complement = 'a';   
    }

    //return the complement
    return complement;
      
}

// returns the complementary DNA string of a valid DNA string
string complementStr(string DNA){
    int length = DNA.length(),i=0;
    string complement;
    //iterates through each character,adding complement letter into the string complement.
    while(i < length){
        complement += complementChar(DNA[i]);
        i++;
    }
    //returns the string complement
    return complement;
    
    
    
    
}

//takes an RNA string and returns a vector of codons contained in an RNA string where each codon is a string of three nucleotides
vector<string> codons(string RNA){
    vector<string> codonsV;
    string codon;
    
    // iterates through every char in the string,when it hits record the first ,second and third char into codon, then take that string and push it into the codonsV vector,then set codon back to an empty string and repeat for the next three characters.
    for(int i=0; i <= RNA.length()-1 ; i += 3){
        codon += RNA[i];
        codon += RNA[i+1];
        codon += RNA[i+2];
        // if the length of the RNA sequence has a remainder of one when divided by 3,the last element only has one char
        if (i == RNA.length()-1 && RNA.length()%3 == 1){
            codon = RNA[i];
            codonsV.push_back(codon);
        // if the length of the RNA sequence has a remainder of 2 when divided by 3,the last element only has two chars.
        }
        else if (i == RNA.length()-2 && RNA.length()%3 == 2){
            codon = "";
            codon += RNA[i];
            codon += RNA[i+1];
            codonsV.push_back(codon);
        }
        else{
            codonsV.push_back(codon);
            codon = "";
        }
    }
    
    
    return codonsV;
    
    
}

//takes a vector of codons (created by the function above), and a three character string as arguments and returns the number of times the string exists in the vector.
int numCodon(vector<string> codonV,string codon){
    int counter=0;
    for(int i =0; i < codonV.size(); i++){
        if(codon == codonV[i]){
            counter++;   
        }
    }
    return counter;
    
}
