#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>

/*------------------------------------------------------------------------------------------------------------------------------------
                        						SECTION 1: DATA STRUCTURE FOR K-MERS 
-------------------------------------------------------------------------------------------------------------------------------------*/
struct kmer_s{
	char kmerChars[64]; 
	
};

struct node_s {
	struct kmer_s data;
	struct node_s *next;
};

struct hashtable_s {
	int size;
	struct node_s **table;	
};


/* Create a new hashtable. */
struct hashtable_s *createHashTable( int size ) {
	
	struct hashtable_s *newHashTable = NULL;
	int i;
	
	/* Allocate the table itself. */
	newHashTable = malloc( sizeof( struct hashtable_s ) );
	
	/* Allocate pointers to the head nodes. */
	newHashTable->table = malloc( sizeof( struct hashtable_s * ) * size );
	
	/* Setting each index containing a head node to Null*/
	for( i = 0; i < size; i++ ) {
		newHashTable->table[i] = NULL;
	}
	
	newHashTable -> size = size;
	
	return newHashTable;
}

/*Create new node to be entered into hashTable*/
struct node_s *createNode(struct kmer_s chars){
	struct node_s *newNode;
	newNode = malloc(sizeof(struct node_s));
	newNode->data=chars;
	newNode->next = NULL;
	return newNode;
}

/*Adding a node into a hash table. Utilizes @createNode*/
void addNode(struct hashtable_s *targetTable,int index,struct kmer_s chars){
	struct node_s * newEntry = NULL;
	struct node_s * current = NULL;
	struct node_s * previous = NULL;
	newEntry = createNode(chars);
	current = targetTable -> table[index];
	
	while( current != NULL){
		previous = current;
		current = current->next;
	}
	if (current == targetTable->table[index]){ //Add at the beginning of a list
		newEntry->next = current;
		targetTable->table[index] = newEntry;
	}else{ //Add at the end of a list
	 previous->next = newEntry;
	}	
}

void removeNode(struct hashtable_s *targetTable,int index,struct kmer_s chars){
	struct node_s * toRemove = NULL;
	struct node_s * current = NULL;
	struct node_s * previous = NULL;
	toRemove = createNode(chars); //create an instance of what we WANT TO REMOVE
	current = targetTable->table[index];
	while(current!=NULL){
		if((strcmp(current->data.kmerChars,toRemove->data.kmerChars)==0)){
			if(previous == NULL){ //At the start of the list
			 targetTable->table[index] = current -> next;
			return;
			}else{ // Any where else in the list
				 previous->next = current -> next;
			return;
			}
		}
	previous = current;
	current = current -> next;
	}
	//printf("No matches found\n");
	return;
}

int findNode(struct hashtable_s *targetTable,int index,struct kmer_s chars){
	struct node_s * toFind = NULL;
	struct node_s * current = NULL;
	toFind = createNode(chars); //create an instance of what we WANT TO REMOVE
	current = targetTable->table[index];
	while(current!=NULL){
		if((strcmp(current->data.kmerChars,toFind->data.kmerChars)==0)){
			return 1;
		}
	
		current = current -> next;
	}
	//printf("No matches found\n");
	return -1;
}

/* Simple method to print the nodes along with how many there are */
void printTable(struct hashtable_s *targetTable, FILE *outFile){
	int n = 0;
	int i;
	for(i=0;i<targetTable->size;i++){
		struct node_s * current = NULL;
		current = targetTable -> table[i];
		while(current != NULL){
			fprintf(outFile,"%s\n", current->data.kmerChars);
			current = current->next;
		}
	}
}

int hashCode(char *str){ //djb2 hashcode
    int hashCode = 5381;
    int c;

    while (c = *str++){
        hashCode = ((hashCode << 4) + hashCode) + c; 
	}
	if(hashCode<0){
		hashCode = hashCode * -1;
	}
    return hashCode;
}

/*------------------------------------------------------------------------------------------------------------------------------------
                        						SECTION 2: DEALING WITH THE FILES 
-------------------------------------------------------------------------------------------------------------------------------------*/


void kmerKram(struct hashtable_s *newTable, FILE *fp,int addOrRemove){
	char c;
	char completeKmer[20];
	int kmerStartFlag = 0;
	int kCount = 0;
	while(1){
		c = fgetc(fp);
		if(c==EOF){
			return;
		}
		if(isalpha(c)){ // The current character is a letter
			if(kmerStartFlag == 0){ // Found the start of a kmer
				kmerStartFlag = 1;
				completeKmer[kCount] = c;
				kCount++;
			}else{ 					// Continue kmer construction
				completeKmer[kCount] = c;
				kCount++;
			}
		}else{ 						// The current character is NOT a letter
			if(kmerStartFlag == 0){ // Ignore character if it's not the end of a kmer
				continue;
			}else{					// Found a character AFTER the end of a kmer
				kmerStartFlag = 0;
				completeKmer[kCount] = '\0';
				kCount = 0;
				struct kmer_s newKmer;
				strcpy(newKmer.kmerChars,completeKmer);
				int index = hashCode(newKmer.kmerChars);
				index = index % newTable->size;
				if(addOrRemove==1){
					addNode(newTable,index,newKmer);
				}else if(addOrRemove==-1){
					removeNode(newTable,index,newKmer);
				}
			}
		}
	}
}

void detectKmer(char *DNAsequence,FILE *outFile,int kmerLength, struct hashtable_s *table,int *uPtr){
	int i;
	int j;
	struct kmer_s testKmer;
	for(i=0;i<strlen(DNAsequence);i++){
		int kCount = 0;
		if(i+kmerLength>strlen(DNAsequence)){
			break;
		}
		for(j=0;j<kmerLength;j++){
			testKmer.kmerChars[kCount]=DNAsequence[j+i];
			kCount++;
		}
		testKmer.kmerChars[kCount] = '\0';
		int index = hashCode(testKmer.kmerChars);
		index = index % table->size;
		if(findNode(table,index,testKmer)==1){
			fprintf(outFile,">%i\n",*uPtr);
			fprintf(outFile,"%s\n", DNAsequence);
			(*uPtr)++;
			break;
		}
	}
}

void skimThroughFasta(FILE *fp,FILE *outFile,int kmerLength, struct hashtable_s *table){
	
	char buffer[512];
	char completeDNA[512];
	int countDNA=0;
	int i;
	int uniqueCounter = 1;
	int *uPtr = &uniqueCounter; 
	char c;
	fgets(buffer,512,fp);
	while(1){
		c = fgetc(fp);
		switch(c){
			case EOF:
				completeDNA[countDNA]='\0';
				detectKmer(completeDNA,outFile,kmerLength,table,uPtr);
				return;
				
			case '>':
				fgets(buffer,512,fp);
				completeDNA[countDNA]='\0';
				detectKmer(completeDNA,outFile,kmerLength,table,uPtr);
				countDNA=0;	
				continue;
			
			case '\r':
				continue;
				
			case '\n':
				continue;
				
			default:
				completeDNA[countDNA]=c;
				countDNA++;
				continue;
		}
	}
}

int fileCounter (FILE *tmpfp){
	int lines = 0;
	char ch;
	while(!feof(tmpfp)){
	ch = fgetc(tmpfp);
		if(ch == '\n'||ch == '\r'){
    		lines++;
 		}
	}
	return lines;
}

void collisionCalculation(struct hashtable_s *testHashTable){
	int collisionArray[20] = { 0 };
	int i;
	for(i=0;i<testHashTable->size;i++){
		int number = 0;
		struct node_s * current = NULL;
		current = testHashTable -> table[i];
		while(current != NULL){
			current = current->next;
			number++;
		}	
		collisionArray[number] +=1;
	}
	for(i=0;i<20;i++){
		printf("Number of indexes with %i: %i \n", i, collisionArray[i]);
	}
}

main(int argc, char *argv[]){
	
	int tabSize = pow(4, atoi(argv[3])) ;
	char hashSize = (char)tabSize;
	pid_t parentpid;
	pid_t pid[1];
	int status;	
	
	char jfCountFile1[256];
	char jfCountFile2[256];
	char jfDumpFile1[256];
	char jfDumpFile2[256];
	char fileName1[256];
	char fileName2[256];
	
	sprintf(jfCountFile1,"time jellyfish count %s -m %s -s %i -o %s-%s-jfCount",argv[1],argv[3],tabSize, argv[1],argv[3]);
	sprintf(jfCountFile2,"time jellyfish count %s -m %s -s %i -o %s-%s-jfCount",argv[2],argv[3],tabSize, argv[2],argv[3]);
	sprintf(jfDumpFile1,"time jellyfish dump %s-%s-jfCount -c -o %s-%s-Kmers",argv[1],argv[3],argv[1],argv[3]);
	sprintf(jfDumpFile2,"time jellyfish dump %s-%s-jfCount -c -o %s-%s-Kmers",argv[2],argv[3],argv[2],argv[3]);
	sprintf(fileName1,"%s-%s-Kmers",argv[1],argv[3]);
	sprintf(fileName2,"%s-%s-Kmers",argv[2],argv[3]);
	printf("%s\n%s\n",jfCountFile1,jfCountFile2);
	printf("%s\n%s\n%s\n%s\n",jfDumpFile1,jfDumpFile2,fileName1,fileName2);
	
	
	if((pid[0] = fork())<0){
		perror("fork");
		abort();
	}else if(pid[0]==0){
		printf("Executing jellyfish functions on %s with %s-mers...\n",argv[1],argv[3]);
		system(jfCountFile1);
		printf("...Count file for %s completed\n",argv[1]);
		system(jfDumpFile1);
		printf("... jellyfish dump completed for %s with %s-mers.\n",argv[1],argv[3]);
		exit(0);
	}
	if((pid[1] = fork())<0){
		perror("fork");
		abort();
	}else if(pid[1]==0){
		printf("Executing jellyfish functions on %s with %s-mers...\n",argv[2],argv[3]);
		system(jfCountFile2);
		printf("...Count file for %s completed\n",argv[2]);
		system(jfDumpFile2);
		printf("... jellyfish dump completed for %s with %s-mers.\n",argv[2],argv[3]);
		exit(0);
	}
	
	while ((parentpid = wait(&status)) > 0);
	puts("Both processes are completed.");
	

	FILE *fp1 = fopen(fileName1,"r");
	FILE *fp2 = fopen(fileName2,"r");
	
	int lineCount = fileCounter(fp2);
	fseek(fp2, 0, SEEK_SET);
	
	struct hashtable_s *uniqueTable = createHashTable(1.5 * lineCount);
	
	kmerKram(uniqueTable,fp2,1);
	kmerKram(uniqueTable,fp1,-1);
	
	char outFileName[256];
	sprintf(outFileName,"%s-VS-%s-%sMERS.fasta",argv[1],argv[2],argv[3]);
	FILE *outFile = fopen(outFileName,"w");
	FILE *sequences = fopen(argv[2],"r");
	fseek(fp2, 0, SEEK_SET);
	
	skimThroughFasta(sequences, outFile, atoi(argv[3]),uniqueTable);
	
}

