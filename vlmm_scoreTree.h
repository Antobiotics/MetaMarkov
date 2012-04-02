#ifndef SEQAN_HEADER_VLMM_SCORETREE_H
#define SEQAN_HEADER_VLMM_SCORETREE_H

namespace SEQAN_NAMESPACE_MAIN
{

// this is the specialisation for the VLMM  used as a score tree
struct ScoreTree{
	};

//create a score tree from two vlmms
template<typename TAlphabet, typename TCargo, typename TVLMMSpec1, typename TVLMMSpec2>
inline void
buildUnionScoreTree(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec1 > > > > & foreground,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec2 > > > > & background)
{
//cout << "Rentre..";
typedef Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec1 > > > > TGraph;
typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef typename Size<TAlphabet>::Type TSize;
TSize alphaSize = ValueSize<TAlphabet>::VALUE;
typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
TVertexDescriptor nilVal =  getNil<TVertexDescriptor>();
//traverse both models (foreground and background) in parallel using their reverse suffix links
TVertexDescriptor foreRoot = getRoot(foreground);
TVertexDescriptor backRoot = getRoot(background);
//should be root
TVertexDescriptor scoreRoot = addAdditionalVertex(scoreTree);
SEQAN_ASSERT(scoreRoot==0)
assignRoot(scoreTree,scoreRoot);
setMarked(scoreTree,scoreRoot,true);
//	cout << "Plante pour GetScoreDistribution"<<endl;
getScoreDistribution(scoreTree,foreground,background,scoreRoot,foreRoot,backRoot);
TVertexDescriptor foreDummy,backDummy;
//cout << "boucle...";
for(TSize i = 0;i<alphaSize;++i){
	foreDummy =	getReverseSuffixLink(foreground,foreRoot,i);
	backDummy =	getReverseSuffixLink(background,backRoot,i);
	
	if(foreDummy != nilVal && backDummy != nilVal && isMarked(foreground,foreDummy) && isMarked(background,backDummy)){
		setScoreRecursively(scoreTree,foreground,background,scoreRoot,foreDummy,backDummy,i);
	}
	if(foreDummy != nilVal && isMarked(foreground,foreDummy) && (backDummy == nilVal || !isMarked(background,backDummy))){
		setScoreForeRecursively(scoreTree,foreground,background,scoreRoot,foreDummy,backRoot,i);
	}
	if( (foreDummy == nilVal || !isMarked(foreground,foreDummy)) && backDummy != nilVal && isMarked(background,backDummy)){
		setScoreBackRecursively(scoreTree,foreground,background,scoreRoot,foreRoot,backDummy,i);
	}

}

//cout << "Fin"<<endl;
}

//create a score tree from two vlmms by recursively going down in both
template<typename TAlphabet, typename TCargo,typename TVLMMSpec1, typename TVLMMSpec2 ,typename TVertexDescriptor,typename TPos>
inline void
setScoreRecursively(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec1 > > > > & foreground,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec2 > > > > & background,
					TVertexDescriptor &scoreNode,
					TVertexDescriptor &foreNode,
					TVertexDescriptor &backNode,
					TPos pos)
{
typedef typename Size<TAlphabet>::Type TType;
TType alphaSize = ValueSize<TAlphabet>::VALUE;
TVertexDescriptor newNode = addAdditionalVertex(scoreTree);
setReverseSuffixLink(scoreTree,scoreNode,newNode,pos);
setSuffixLink(scoreTree,newNode,scoreNode);
setMarked(scoreTree,newNode,true);
getScoreDistribution(scoreTree,foreground,background,newNode,foreNode,backNode);
TVertexDescriptor foreDummy,backDummy;
TVertexDescriptor nilVal =  getNil<TVertexDescriptor>();
for(TType i = 0;i<alphaSize;++i){
	foreDummy =	getReverseSuffixLink(foreground,foreNode,i);
	backDummy =	getReverseSuffixLink(background,backNode,i);
	if(foreDummy != nilVal && backDummy != nilVal && isMarked(foreground,foreDummy) && isMarked(background,backDummy)){
		setScoreRecursively(scoreTree,foreground,background,newNode,foreDummy,backDummy,i);
	}
	if(foreDummy != nilVal && isMarked(foreground,foreDummy) && (backDummy == nilVal || !isMarked(background,backDummy))){
		setScoreForeRecursively(scoreTree,foreground,background,newNode,foreDummy,backNode,i);
	}
	if( (foreDummy == nilVal || !isMarked(foreground,foreDummy)) && backDummy != nilVal && isMarked(background,backDummy)){
		setScoreBackRecursively(scoreTree,foreground,background,newNode,foreNode,backDummy,i);
	}

}
}

//create a score tree from two vlmms by recursively going down in the foreground model
template<typename TAlphabet, typename TCargo, typename TVLMMSpec1, typename TVLMMSpec2 ,typename TVertexDescriptor,typename TPos>
inline void
setScoreForeRecursively(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec1 > > > > & foreground,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec2 > > > > & background,
					TVertexDescriptor &scoreNode,
					TVertexDescriptor &foreNode,
					TVertexDescriptor &backNode,
					TPos pos)
{
typedef typename Size<TAlphabet>::Type TType;
TType alphaSize = ValueSize<TAlphabet>::VALUE;
TVertexDescriptor newNode = addAdditionalVertex(scoreTree);
setReverseSuffixLink(scoreTree,scoreNode,newNode,pos);
setSuffixLink(scoreTree,newNode,scoreNode);
setMarked(scoreTree,newNode,true);
getScoreDistribution(scoreTree,foreground,background,newNode,foreNode,backNode);
TVertexDescriptor foreDummy;
TVertexDescriptor nilVal =  getNil<TVertexDescriptor>();
for(TType i = 0;i<alphaSize;++i){
	foreDummy =	getReverseSuffixLink(foreground,foreNode,i);
	if(foreDummy != nilVal && isMarked(foreground,foreDummy))
			setScoreForeRecursively(scoreTree,foreground,background,newNode,foreDummy,backNode,i);
	

}
}

//create a score tree from two vlmms by recursively going down in the background model
template<typename TAlphabet, typename TCargo, typename TVLMMSpec1, typename TVLMMSpec2 ,typename TVertexDescriptor,typename TPos>
inline void
setScoreBackRecursively(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec1 > > > > & foreground,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec2 > > > > & background,
					TVertexDescriptor &scoreNode,
					TVertexDescriptor &foreNode,
					TVertexDescriptor &backNode,
					TPos pos)
{
typedef typename Size<TAlphabet>::Type TType;
TType alphaSize = ValueSize<TAlphabet>::VALUE;
TVertexDescriptor newNode = addAdditionalVertex(scoreTree);
TVertexDescriptor nilVal =  getNil<TVertexDescriptor>();
setReverseSuffixLink(scoreTree,scoreNode,newNode,pos);
setSuffixLink(scoreTree,newNode,scoreNode);
setMarked(scoreTree,newNode,true);
getScoreDistribution(scoreTree,foreground,background,newNode,foreNode,backNode);
TVertexDescriptor backDummy;
for(TType i = 0;i<alphaSize;++i){
	backDummy =	getReverseSuffixLink(background,backNode,i);
	if( backDummy != nilVal && isMarked(background,backDummy) ){
		setScoreBackRecursively(scoreTree,foreground,background,newNode,foreNode,backDummy,i);
	}

}
}

// Node access on the score tree

//create a score tree from two vlmms by recursively going down in both
template<typename TAlphabet, typename TCargo, typename TVLMMSpec1, typename TVLMMSpec2 ,typename TVertexDescriptor>
inline void
getScoreDistribution(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec1 > > > > & foreground,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec2 > > > > & background,
					TVertexDescriptor &scoreNode,
					TVertexDescriptor &foreNode,
					TVertexDescriptor &backNode)
{
	//cout << "GetScore IN..";
typedef typename Size<TAlphabet>::Type TType;
TType alphaSize = ValueSize<TAlphabet>::VALUE;
for(TType i = 0;i<alphaSize;++i){
	//cout<< "F/B: "<<getProbability(foreground,foreNode,i)<<" /" <<getProbability(background,backNode,i);
	setScore(scoreTree,scoreNode,i, rounded(log(getProbability(foreground,foreNode,i)/getProbability(background,backNode,i))) );
//	cout <<" Result: " <<getScore(scoreTree,scoreNode,i)<<endl;
}
//	cout << "OUT !"<<endl;
}

double rounded(double number){
	return number;
	//return (double)floor(number+.5);

}

template<typename TCargo,typename TAlphabet,typename TChar ,typename TVertexDescriptor>
inline double
getScore(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &vlmm,
			   TVertexDescriptor &father,
			   TChar pos)
{
	TAlphabet letter(pos);
	return value(vlmm.data_probability_vector[father],(int)letter);
}

template<typename TCargo,typename TAlphabet,typename TChar, typename TVertexDescriptor>
inline void
setScore(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &vlmm,
			   TVertexDescriptor &father,
			   TChar pos,
			   double score)
{
	TAlphabet letter(pos);
	value(vlmm.data_probability_vector[father],(int)letter) = score;
}


/*
* Score retrieval
*/
// returns the summed scores for a text, this is called normalized score Mode (ns)
// in case that the score tree is build via likelihood ratios
// this is exactly the log-likelihood ratio (log(H1(s)/H0(s))) of the score tree for the whole sequence  
template<typename TAlphabet,typename TCargo>
inline double
estimateScore( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &vlmm,
					String<TAlphabet> &text)
{
double result = 0;
	typedef typename Iterator<String<TAlphabet>, Rooted  >::Type TIter;
	TIter it = begin(text,Rooted());

	for(;!atEnd(it);goNext(it))
	{
		result += getProbabilityForLongestContext(vlmm,it);
		//std::cout <<" prob for letter: "<<value(it)<< " is: "<<getProbabilityForLongestContext(vlmm,it)<<endl;
	}
return result;
}

	
template<typename TAlphabet,typename TCargo>
inline vector<double>
estimateScore2( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &vlmm,
				  String<TAlphabet> &text)
	{
		vector<double> ResVect;
		double result = 0;
		typedef typename Iterator<String<TAlphabet>, Rooted  >::Type TIter;
		TIter it = begin(text,Rooted());
		
		for(;!atEnd(it);goNext(it))
			//log ou pas log ?
			//getProbabilityForLongestContext ou pas?
		{	double tt=getProbabilityForLongestContext(vlmm,it);
			result += tt;
			if(result>0){
				
			
			ResVect.push_back(result);
			}
			else{
				ResVect.push_back(0);
			}
			//std::cout <<" prob for letter: "<<value(it)<< " is: "<<getProbabilityForLongestContext(vlmm,it)<<endl;
		}
		return ResVect;
	}
	
	
	
	
	
// gets the scores for a text, this is called single base score Mode (sbs)
template<typename TAlphabet,typename TCargo>
inline void
estimateScoreValues( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &vlmm,
					String<TAlphabet> &text,
					String<double> &values)
{
	resize(values,length(text));
	typedef typename Iterator<String<TAlphabet>, Rooted>::Type TIter;
	typedef typename Iterator<String<double>, Standard >::Type TIterDouble;
	TIter it = begin(text,Rooted());
	TIterDouble itVal = begin(values,Standard());
	for(;!atEnd(it);goNext(it), goFurther(itVal,1))
		assignValue(itVal,getProbabilityForLongestContext(vlmm,it));

	return;
}


// gets the scores for a text and writes directly into a file, in single base score Mode (sbs)
template<typename TAlphabet,typename TCargo,typename TFile>
inline void
estimateScoreValues( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &vlmm,
					String<TAlphabet> &text,
					TFile &outFile)
{
	typedef typename Iterator<String<TAlphabet>, Rooted  >::Type TIter;
	TIter it = begin(text,Rooted());
	for(;!atEnd(it);goNext(it)){
			_streamPut(outFile,'\t');
			_streamPutDouble(outFile,getProbabilityForLongestContext(vlmm,it));
	}
	return;
}



template<typename TAlphabet,typename TCargo>
inline void
getBestScores( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > const &vlmm,
					String<double> &bestScores,
					String<String<TAlphabet> > &scoreLabels)
{

	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<ScoreTree> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef pair<double,TVertexDescriptor> entry;
	std::priority_queue<entry> pQueue;
	//std::priority_queue<double> pQueue;
	//typedef Graph<Directed<> > TGraphD;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	//typedef typename Position<String<AutomatonEdgeArray<TEdge, TAlphabet> > >::Type TPos;
	//typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> >, Rooted >::Type TIterConst;
	//TIterConst it = begin(vlmm.data_vertex);
	TVertexIterator it(vlmm);
	for(goBegin(it);!atEnd(it);goNext(it)) {
		//TPos pos = position(it);
		//if (!idInUse(vlmm.data_id_managerV, position(it))) continue;
		TVertexDescriptor sourceVertex =*(it);//position(it);

		for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i){
			pQueue.push(entry(value(vlmm.data_probability_vector[sourceVertex],i),sourceVertex));
		}

	}
	//size of bestScores is the number of Scores to return
	for(unsigned j = 0;j<length(bestScores);++j){
		entry dummy = pQueue.top();
		bestScores[j] = dummy.first;
		TVertexDescriptor father,child=dummy.second;
		for(int i = 0;i<ValueSize<TAlphabet>::VALUE;++i){
			double score = value(vlmm.data_probability_vector[child],i);
			if(dummy.first == score ){
					appendValue(scoreLabels[j],(TAlphabet)i);
					while(!isRoot(vlmm,child)){
						father = getSuffixLink(vlmm,child);
						appendValue(scoreLabels[j],getReverseSuffixLinkCharacter(vlmm,father,child));
						child = father;
					}

				break;
			}
		}
		pQueue.pop();
	}
}

// create the expected score of a vlmm scoring tree under a vlmm background model
// we assume that the score tree has as its topology the union of the foreground and background model
// so it is the restrictive component, i.e. no branch in the background model can be longer than in the score Tree
// where the opposite might be the case
template<typename TAlphabet, typename TCargo, typename TVLMMSpec>
inline double
expectedValue(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
			  Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > & background){
typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TVLMMSpec> > > > TGraph;
typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
String<unsigned> longestLeaves;
resize(longestLeaves,1);
String<String<TAlphabet> > leafLabels;
resize(leafLabels,1);
getLongestLeaves(scoreTree,longestLeaves,leafLabels);
// the longest branch in the background model
// means we calculate the expected value for all words of length maxDepth+1
unsigned maxDepth = length(leafLabels[0]);
std::cout<<" the longest leaf in the score tree as length:"<<maxDepth<<endl;
TVertexDescriptor scoreRoot=getRoot(scoreTree),backRoot=getRoot(background);
double E = getExpectedValueRecursively(scoreTree,background,scoreRoot,backRoot,maxDepth);
cout << "The expected value E under the given probability model is E:"<<E<<endl<<endl;
return E;
}


//create a score tree from two vlmms by recursively going down in both
template<typename TAlphabet, typename TCargo, typename TVLMMSpec ,typename TVertexDescriptor>
inline double
getExpectedValueRecursively(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > & background,
					TVertexDescriptor &scoreNode,
					TVertexDescriptor &backNode,
					unsigned restDepth)
{
typedef typename Size<TAlphabet>::Type TType;
TType alphaSize = ValueSize<TAlphabet>::VALUE;
TVertexDescriptor scoreDummy,backDummy;
TVertexDescriptor nilVal =  getNil<TVertexDescriptor>();
double result=0,howManyChild = pow((double)alphaSize,(double)restDepth);
for(TType i = 0;i<alphaSize;++i){
	scoreDummy =	getReverseSuffixLink(scoreTree,scoreNode,i);
	backDummy =	getReverseSuffixLink(background,backNode,i);
	if(scoreDummy != nilVal && backDummy != nilVal && isMarked(scoreTree,scoreDummy) && isMarked(background,backDummy)){
		result += getExpectedValueRecursively(scoreTree,background,scoreDummy,backDummy,restDepth-1);
	}
	if(scoreDummy != nilVal && isMarked(scoreTree,scoreDummy) && (backDummy == nilVal || !isMarked(background,backDummy))){
		result += getExpectedValueScoreTreeRecursively(scoreTree,background,scoreDummy,backNode,restDepth-1);
	}
	else{ // get the probabilities for that letter
		result += howManyChild*getProbability(background,backNode,i)*getScore(scoreTree,scoreNode,i);
	}


}
return result;
}

//create a score tree from two vlmms by recursively going down in the foreground model
template<typename TAlphabet, typename TCargo, typename TVLMMSpec ,typename TVertexDescriptor>
inline double
getExpectedValueScoreTreeRecursively(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
					Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > & background,
					TVertexDescriptor &scoreNode,
					TVertexDescriptor &backNode,
					unsigned restDepth)
{
typedef typename Size<TAlphabet>::Type TType;
TType alphaSize = ValueSize<TAlphabet>::VALUE;
double result=0,howManyChild = pow((double)alphaSize,(double)restDepth);
TVertexDescriptor scoreDummy;
TVertexDescriptor nilVal =  getNil<TVertexDescriptor>();
for(TType i = 0;i<alphaSize;++i){
	scoreDummy =	getReverseSuffixLink(scoreTree,scoreNode,i);
	if(scoreDummy != nilVal && isMarked(scoreTree,scoreDummy))
			result += getExpectedValueScoreTreeRecursively(scoreTree,background,scoreDummy,backNode,restDepth-1);
	else{ // get the probabilities for that letter
		result += howManyChild*getProbability(background,backNode,i)*getScore(scoreTree,scoreNode,i);
	}
}
return result;
}


template<typename TAlphabet,typename TCargo>
inline void
getLeafScores( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > >  &vlmm,
					String<double> &Leaves)
{

	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<ScoreTree> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TAlph;
	TAlph alphaSize = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor dummy,root=getRoot(vlmm),nilVal =  getNil<TVertexDescriptor>();
	typedef pair<unsigned,TVertexDescriptor> entry;
	std::priority_queue<entry> pQueue;
	unsigned count = 0,depth=0;
	for(TAlph i = 0;i<alphaSize;++i){
		dummy =	getReverseSuffixLink(vlmm,root,i);
		if( dummy != nilVal && isMarked(vlmm,dummy) ){
			traverseReverseSuffixLinks(vlmm,dummy,pQueue,depth+1);
			++count;
		}
	}
	if(count==0){
		// only rhe root node will be added
		entry set(0,root);
		pQueue.push(set);
	}
	while(!pQueue.empty()){
		entry dummy = pQueue.top();
		for(TAlph i = 0;i<alphaSize;++i){
			double score = getScore(vlmm,dummy.second,i);
			append(Leaves,score);
		}
		pQueue.pop();
	}
}

/*************
*
* Local Score
*
**************/
	
template<typename TAlphabet, typename TCargo>
double
maximumScoreMulti(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
		unsigned maxDepth,
		String<TAlphabet> &sequence,
		pair<unsigned,unsigned> &stats)
	{

		
		double largeNegativeNumber=-100000000;
		double optimal=largeNegativeNumber;
		unsigned L;
		unsigned i=1;
		//Je n'ai aucune idee de valeur pour T.
		double T=10;
		
		String<double> max_before,max_so_far,max_here;
		String<unsigned> start_before;
		resize(max_before,length(sequence)+1);
		resize(start_before,length(sequence)+1);
		resize(max_here,length(sequence)+1);
		resize(max_so_far,length(sequence)+1);
		
		max_before[0]=largeNegativeNumber;
		max_here[0]=largeNegativeNumber;
		max_so_far[0]=largeNegativeNumber;
		//Bidouille:
		start_before[0]=1;
		
		String<double> scores,subScore;
		resize(subScore,maxDepth+1);
		resize(scores,maxDepth+1);
		
		typedef typename Iterator<String<TAlphabet>,Rooted>::Type TIter;
		TIter it = begin(sequence,Rooted());
		cout << "--------------------------"<<endl;
		for(;!atEnd(it);goNext(it))
		{	
			if (i>1) {
				
			
			
			L=getScoreForAllContext(scoreTree,it,scores);
			
			max_here[i]=largeNegativeNumber;
			///cout <<"Partie 1:" <<endl;
			unsigned width=-1;
			if ((i+1)<L) {
				width=i+1;
			}
			else{
				width=L;
			}
			//cout << "width: "<< width<< endl;
			for(unsigned j=1;j<=width;++j)
			{
				subScore[j] = scores[j] + subScore[j-1];
				//cout << subScore[j]<<endl;
				if(optimal < subScore[j]){
					optimal = subScore[j];
					
					
				}
				if(optimal*max_so_far[i-j-1]>=max_here[i]){
					max_here[i]=optimal*max_so_far[i-j-1];
					start_before[i]=j;
			//		cout << "j: "<< j <<endl;
				}
				//Pas sur :
				else{
					start_before[i]=i;
			//		cout<< "passe par la"<<endl;
				}
				
				
			}
			//cout<<"Partie 2:"<<endl;
			if ((i>=L)) {
				
				max_before[i]=max_before[i-1]+scores[L];// *probs[L-1] ou + autre chose
				
			//	cout << "la ?"<<endl;
				
			//	cout << start_before[i]<<endl;
				
				double temp=max_before[i]*max_so_far[start_before[i]-1];

				if (temp>=max_here[i]) {
					max_here[i]=temp;
				}
				
			}
			//cout<<"Partie 3:"<<endl;
			if (max_so_far[i-1]>=max_here[i]-T) {
				max_so_far[i]=max_so_far[i-1];
			}
			else {
				max_so_far[i]=max_here[i]-T;
			}
			//cout<<"Partie 4:"<<endl;
			if (i>=(L-1)) {
				if (subScore[L-1]>=max_before[i]) {
					max_before[i]=subScore[L-1];
				}
			}
			}
			i++;
		}
		/*
		cout << "max_before"<<endl;
		for (int i=0; i<length(max_before); i++) {
			cout << i <<":" <<max_before[i]<< " ";
		}
		cout << endl;
		cout << "max_here"<<endl;
		for (int i=0; i<length(max_here); i++) {
			cout << i <<":" <<max_here[i]<< " ";
		}
		cout << endl;
		cout << "max_so_far"<<endl;
		for (int i=0; i<length(max_so_far); i++) {
			cout << i <<":" <<max_so_far[i]<< " ";
		}
		cout << endl;
		
		for (int i=0; i<length(start_before); i++) {
			cout << i <<":" <<start_before[i]<< " ";
		}
		cout << endl;*/
		cout << max_so_far[i]<<endl;
		return(max_so_far[i]);
	}
//look for the maximum scoring segment among all possible segments of the input sequence
//where the scoring distribution comes from a vlmm scoring tree
template<typename TAlphabet, typename TCargo>
double
maximumScore(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > & scoreTree,
			 unsigned maxDepth,
			 String<TAlphabet> &sequence,
			 pair<unsigned,unsigned> &stats)
{
	//cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	unsigned L,j,step=1,length1=0,prevLength=0;
	double optimal=-1000,previous=0;//,formeropt;
	String<double> scores,subScore; // the scores are listed with increasing Context
	resize(scores,maxDepth+1);		// i.e. subScore(0) = P(s_i) or subScore(1) = P(s_i|s_i-1)
	resize(subScore,maxDepth+1);
	typedef typename Iterator<String<TAlphabet>,Rooted >::Type TIter;
	TIter it = begin(sequence,Rooted());
	unsigned index = 0;
	for(;!atEnd(it);goNext(it))
	{
		// remember the length of longest leaf in L and all the scores up to length L in array scores
		L = getScoreForAllContext(scoreTree,it,scores);
		cout <<"L:"<<L<<" ";
	/*	for(unsigned h=0;h<=length(scores);++h)
			cout << "sco["<<h<<"]: "<<scores[h]<<" ";
		cout << endl;
		// now process all prefix scores of the current string position s_i and check if any is better than the optimal yet
		cout << length(scores)<< " " << L << endl;*/
		for(j=1;j<=L;++j)
		{
			subScore[j] = scores[j] + subScore[j-1];
			cout << subScore[j]<<endl;
			if(optimal < subScore[j]){
				optimal = subScore[j];
				index =step;
				length1 = j+1;
			}
			/*optimal = max(optimal,subScore[j]);
			if(formeropt< optimal){
				formeropt = optimal;
				index = step;
			}*/
		}
		subScore[0] = scores[0];

		if(optimal < subScore[0]){
				optimal = subScore[0];
				index =step;
				length1 = 1;
		}
		/*optimal = max(optimal,subScore[0]);
		if(formeropt< optimal){
				formeropt = optimal;
				index = step;
			}*/
		//for(unsigned h=0;h<=L;++h)
		//	cout << "sub["<<h<<"]: "<<subScore[h]<<" ";
		// always add the actual longest score to the extenion of the previous best score
		previous = previous + scores[L] ;
		++prevLength;
		if(optimal < previous){
				optimal = previous;
				index =step;
				length1 = prevLength;
		}
		//cout << index << " " << length<<endl;
		/*optimal = max(optimal,previous);
			if(formeropt< optimal){
				formeropt = optimal;
				index = step;
			}*/
		// check the best score. Either the best score starts at the longest subscore[L] or is a continuation
		// of a segment longer than L
		if(previous < subScore[L]){
			previous = subScore[L];
			prevLength =L+1;
		}
		//previous = max(subScore[L],previous);
		++step;
	//	cout<< "o: "<<optimal<<" p: "<<previous<<endl;
	}
	//cout <<"the optimal sequence was found ending at position "<<index<<" with length: "<<length<<endl;
	// first = end position of local score segment
	stats.first = index;
	// second is length of local score segment
	stats.second = length1;
	//cout <<"opt: " << optimal<<endl;
	//cout <<"prev: " << previous<<endl;
	if(optimal >0)
		return optimal;
	else{
		stats.first = 0;
		stats.second=0;
		return 0;
	}
	
		
}



template<typename TAlphabet,typename TCargo,typename TVLMMSpec,typename TIter>
inline unsigned
getScoreForAllContext( Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < TVLMMSpec > > > > &vlmm,
					TIter &it ,
					String<double> &scores)
{
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<VLMM<TVLMMSpec> > > > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TIter copy = it;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor node = getRoot(vlmm);
	unsigned depth = 0;
	scores[0] = getProbability(vlmm,node,value(it));
	while(!atBegin(copy) ){
		goPrevious(copy);
	

		if(getReverseSuffixLink(vlmm,node,value(copy)) != nilVal)
		{
				node = getReverseSuffixLink(vlmm,node,value(copy));

				if(isMarked(vlmm,node)){
						++depth;
						//cout <<"letter:" <<value(it)<<"  getProb: "<<getProbability(vlmm,node,value(it))<<" log(getProb): "<<log(getProbability(vlmm,node,value(it)))<<endl;
						scores[depth] = getProbability(vlmm,node,value(it));
				}
				else
						break;
		}
		else
			break;
		
	}
	//return getProbability(vlmm,node,value(it));
	return depth;
}

/****
* save the score Tree
*
*****/
template<typename TFile, typename TCargo >
inline void
writeHead(Graph<Automaton<Dna, TCargo , WordGraph < VLMM < ScoreTree > > > > &vlmm,
	   TFile & target)
{	
	String<unsigned> longestLeaves;
	resize(longestLeaves,1);
	String<String<Dna> > leafLabels;
	resize(leafLabels,1);
	getLongestLeaves( vlmm,longestLeaves,leafLabels);
		_streamWrite(target,"VLMM\tScoreTree\t");
		_streamPutInt(target,longestLeaves[0]);
		_streamWrite(target,"\tDna\t");
		_streamPutInt(target,length(vlmm.data_marked));
}

template<typename TFile, typename TCargo >
inline void
writeHead(Graph<Automaton<AminoAcid, TCargo , WordGraph < VLMM < ScoreTree > > > > &vlmm,
	   TFile & target)
{	
	String<unsigned> longestLeaves;
	resize(longestLeaves,1);
	String<String<AminoAcid> > leafLabels;
	resize(leafLabels,1);
	getLongestLeaves( vlmm,longestLeaves,leafLabels);
	_streamWrite(target,"VLMM\tScoreTree\t");
		_streamPutInt(target,longestLeaves[0]);
		_streamWrite(target,"\tAminoAcid\t");
		_streamPutInt(target,length(vlmm.data_marked));
}


/****************************
* read scoring tree
* and score sequences in file
*****************************/

template<typename TInputFile,typename TFile>
inline void
readForLocalScoreEstimation(TInputFile	& file,
						  String<char>	& sequenceFile,
						  TFile			& outFile)
{
	unsigned maxDepth;
	String<char> entry;
	_scanNextEntry(file,entry);
	_scanNextEntry(file,entry);
	if(entry!= "ScoreTree"){
		cerr<<"No ScoreTree can be found in the input file"<<file<<". Maybe incorrect input file format.\n";
		std::exit(1);
	}
	maxDepth=_scanNextIntEntry(file);
	_scanNextEntry(file,entry);

	if(entry == "Dna"){
		Graph<Automaton<Dna, String<Dna> , WordGraph < VLMM < ScoreTree > > > > vlmmDna;
		readGraph(vlmmDna,file);
		estimateLocalScoreOnFile(vlmmDna,sequenceFile,outFile,maxDepth);
	}
	if(entry == "AminoAcid"){
		Graph<Automaton<AminoAcid, String<AminoAcid> , WordGraph < VLMM < ScoreTree > > > > vlmmProtein;
		readGraph(vlmmProtein,file);
		estimateLocalScoreOnFile(vlmmProtein,sequenceFile,outFile,maxDepth);
	}
	else{
		cerr<<"Alphabet in file "<<file<<" not supported for vlmm. Maybe incorrect input file format.\n";
		std::exit(1);
	}
}
//Local score avec max.
// estimates the likelihood for every sequences in the input file and save it in the outout file
template<typename TAlphabet,typename TCargo, typename TFile>
inline void estimateLocalScoreOnSequences(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &scoreTree,
									 String<String<TAlphabet> > &sequences,
									 String<String<char> > &ids,
									 TFile &outFile,
									 unsigned maxDepth)
{
	cout << "!!!local score !!!!"<< endl;
		typedef typename Iterator<String<String<TAlphabet> >, Rooted  >::Type TIter;
		typedef typename Iterator<String<String<char> >, Rooted  >::Type TIterChar;
		TIter it = begin(sequences,Rooted());
		TIterChar id = begin(ids,Rooted());
		pair<unsigned,unsigned> stats;
		double bestScore;
		//put first line for file
		_streamWrite(outFile,"Id\tLocalScore\tEndPosition\tLength\n");
		for(;!atEnd(it);goNext(it),goNext(id))
		{
			bestScore = maximumScore(scoreTree,maxDepth,value(it),stats);
			_streamWrite(outFile,value(id));
			_streamPut(outFile,'\t');
			_streamPutDouble(outFile,bestScore);
			_streamPut(outFile,'\t');
			_streamPutInt(outFile,stats.first);
			_streamPut(outFile,'\t');
			_streamPutInt(outFile,stats.second);
			_streamPut(outFile,'\n');
		}
}

// estimates the likelihood for every sequences in the input file and save it in the out file
template<typename TAlphabet,typename TCargo>
inline void estimateLocalScoreOnSequences(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &scoreTree,
										  std::priority_queue<pair<double,unsigned> > &temp,
										  String<String<TAlphabet> > &sequences,
										  unsigned maxDepth)
{
		typedef typename Iterator<String<String<TAlphabet> >, Rooted >::Type TIter;
		TIter it = begin(sequences,Rooted());
		unsigned count = 0;
		pair<unsigned , unsigned> dummy;

	double ss=0;
	std::priority_queue<pair<double,unsigned> > pQueue;
	vector<double> scores;
		for(;!atEnd(it);goNext(it),++count)
		{
			//double localScore = maximumScore(scoreTree,maxDepth,value(it),dummy);
			double localScore = maximumScoreMulti(scoreTree,maxDepth,value(it),dummy);
			scores.push_back(localScore);
			ss=ss+(localScore);
			
			pair<double,unsigned> entry(localScore,count);
			pQueue.push(entry);
			//cout << count<<endl;

		}
	cout << "fin"<<endl;
	double smax = 0;
	double smin = 0;
	vector<double>::size_type siz = scores.size();
	for (int i= 0; i<siz; i++) {
		if(scores[i]>smax){
			smax=scores[i];
		}
		if (scores[i]<smin) {
			smin=scores[i];
		}
	}
	ss=ss/count;
	unsigned i =0;
	double mid=0;
	double var=0;
	double Amin=-10;
	double Amax=10;
	//fstream fi;
	//fi.open("/Users/gregoirelejay/Desktop/ttt.txt",ios_base::out);
	//fi << "vals"<<endl;
	while(!pQueue.empty()){
		
		double tt=pQueue.top().first;
		
		double shif = Amin+ ((Amax-Amin)/(smax-smin))*(tt-smin);
		
		//fi <<shif<<endl;
		pair<double,unsigned> entry(shif,count-(i+1));
	
		temp.push(entry);
		pQueue.pop();
		mid+=shif;
		var+=shif*shif;
		++i;
	}
	var=var/count;
	mid=mid/count;
	
	cout<< "Variance: "  <<var<<endl;
	cout<< "Esperance: " <<mid<<endl;
	
		
}

// estimates the likelihood for every sequences in the input file and save it in the outout file
template<typename TAlphabet,typename TCargo, typename TFile>
inline void estimateNormalizedScoreOnSequences(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &scoreTree,
									 String<String<TAlphabet> > &sequences,
									 String<String<char> > &ids,
									 TFile &outFile)
{

		typedef typename Iterator<String<String<TAlphabet> >, Rooted  >::Type TIter;
		typedef typename Iterator<String<String<char> >, Rooted>::Type TIterChar;
		TIter it = begin(sequences,Rooted());
		TIterChar id = begin(ids,Rooted());
		double Score;
		//put first line for file
		_streamWrite(outFile,"Id\tNormalizedScore\n");
		for(;!atEnd(it);goNext(it),goNext(id))
		{
			Score = estimateScore(scoreTree,value(it))/length(value(it));
			_streamWrite(outFile,value(id));
			_streamPut(outFile,'\t');
			_streamPutDouble(outFile,Score);
			_streamPut(outFile,'\n');
		}
}


// estimates the likelihood for every sequences in the input file and save it in the out file
template<typename TAlphabet,typename TCargo>
inline void estimateNormalizedScoreOnSequences(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &scoreTree,
										  std::priority_queue<pair<double,unsigned> > &pQueue,
										  String<String<TAlphabet> > &sequences)
{
		typedef typename Iterator<String<String<TAlphabet> >, Rooted  >::Type TIter;
		TIter it = begin(sequences,Rooted());
		unsigned count = 0;
		for(;!atEnd(it);goNext(it),++count)
		{
			double Score = estimateScore(scoreTree,value(it))/length(value(it));
			// put the local score into the priority Queue
			pair<double,unsigned> entry(Score,count);
			pQueue.push(entry);
		}
}

// estimates the likelihood for every sequences in the input file and save it in the outout file
template<typename TAlphabet,typename TCargo, typename TFile>
inline void estimateSingleBaseScoreOnSequences(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &scoreTree,
									 String<String<TAlphabet> > &sequences,
									 String<String<char> > &ids,
									 TFile &outFile)
{

		typedef typename Iterator<String<String<TAlphabet> >,Rooted >::Type TIter;
		typedef typename Iterator<String<String<char> >, Rooted >::Type TIterChar;
		TIter it = begin(sequences,Rooted());
		TIterChar id = begin(ids,Rooted());
		//put first line for file
		for(;!atEnd(it);goNext(it),goNext(id))
		{   
			_streamWrite(outFile,value(id));
			estimateScoreValues( scoreTree,value(it),outFile);
			_streamPut(outFile,'\n');
		}
}
template<typename TAlphabet,typename TCargo>
inline void outputScoreTree(Graph<Automaton<TAlphabet, TCargo , WordGraph < VLMM < ScoreTree > > > > &scoreTree,
									 String<String<TAlphabet> > &sequences,
									 String<String<char> > &ids,
									 String<char> &output,
									 bool sbs,
									 bool ls,
									 bool ns,
									 unsigned maxDepth)
{
fstream outFile;
if(ns){
			openAndAppend(output,".ns",outFile);
			estimateNormalizedScoreOnSequences(scoreTree,sequences,ids,outFile);
			outFile.close();
			std::cout << "wrote normalized score for sequences into: "<<output<<".ns"<<endl;
		}
if(ls){
			openAndAppend(output,".ls",outFile);
			estimateLocalScoreOnSequences(scoreTree,sequences,ids,outFile,maxDepth);
			outFile.close();
			std::cout << "wrote local score for sequences into: "<<output<<".ls"<<endl;
		}
if(sbs){
			openAndAppend(output,".sbs",outFile);
			estimateSingleBaseScoreOnSequences(scoreTree,sequences,ids,outFile);
			outFile.close();
			std::cout << "wrote single base score for sequences into: "<<output<<".sbs"<<endl;
}


}
									  



template<typename TInputFile,typename TFile>
inline void
readVlmmForScoreTree(TInputFile			& familyFile,
					 TInputFile			& backgroundFile,
					 TFile			    & outFile)
{
	// temporary hard coded. File input should be better!!!
	String<char> entry1,entry2;
	_scanNextEntry(familyFile,entry1);
	_scanNextEntry(familyFile,entry1);
	_scanNextEntry(familyFile,entry1);

	_scanNextEntry(backgroundFile,entry2);
	_scanNextEntry(backgroundFile,entry2);
	_scanNextEntry(backgroundFile,entry2);
	
	if(entry1 == "Dna" && entry2 == "Dna"){
		Graph<Automaton<Dna, String<Dna> , WordGraph < VLMM < ContextTree > > > > familyVlmm,backgroundVlmm;
		readGraph(familyVlmm,familyFile);
		readGraph(backgroundVlmm,backgroundFile);
		Graph<Automaton<Dna, String<Dna> , WordGraph < VLMM < ScoreTree > > > > scoreTree;
		buildUnionScoreTree(scoreTree,familyVlmm,backgroundVlmm);
		save(scoreTree,outFile);

	}else
	
		if(entry1 == "AminoAcid" && entry2 == "AminoAcid"){
			Graph<Automaton<AminoAcid, String<AminoAcid> , WordGraph < VLMM < ContextTree > > > > familyVlmm,backgroundVlmm;
			readGraph(familyVlmm,familyFile);
			readGraph(backgroundVlmm,backgroundFile);
			Graph<Automaton<AminoAcid, String<AminoAcid> , WordGraph < VLMM < ScoreTree > > > > scoreTree;
			buildUnionScoreTree(scoreTree,familyVlmm,backgroundVlmm);
			save(scoreTree,outFile);
		}
		else{
		cerr<<"Alphabet in one of the two input files:"<<familyFile<<" or "<<backgroundFile<<" does either not match or maybe inconsistent. Maybe incorrect input file format.\n";
		exit(1);
		}


}

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
