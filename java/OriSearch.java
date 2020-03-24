public class OriSearch {
	public  char DigtoChar(double d){
		char a=' ';
		if(d==0)a='A';
		if(d==1)a='C';
		if(d==2)a='G';
		if(d==3)a='T';
		return a;
	}
	public  String NumbertoString(double a){
		
		String str = "";
		while(a>=4){
			str=DigtoChar(a%4) + str;
			a=(a-(a%4))/4;
			
		}
		return str;
	}
	public  int TextPattern(String text,String pattern){//gives count of pattern repetition
		int i=0,count=0;
		while(i<text.length()-pattern.length()){
            String test = text.substring(i,i+pattern.length());
			if(test.equals(pattern))count=count + 1;
			i++;
			}
	    return count;
	    }
	public String[] FrequentPatterns(String text,int k,int t){ //finds frequent patterns of length k
		
		int max=0;
		int a = 1;int b=1;
		for(int j = 1;j<=k;j++){
			b = b*4;
		}
		int i=b;
		for(int j = 1;j<=k+1;j++){
			a = a*4;
		}
		int[] count = new int[b];
		String[] frequents = new String[1000];
		while(i<b+b){
            String test = NumbertoString(i);
            String testc= ComplementReverse(NumbertoString(i));//System.out.print( i + test+" ");
			count[i-b]=ApproxMatchCount(test,text,t)+ApproxMatchCount(testc,text,t);
			if(count[i-b]>max)max=count[i-b];
			i++;
			}
		//System.out.print(max);
		i=0;int m = 0;
		while(i<count.length){
			    if(count[i]==max){
	            frequents[m]=NumbertoString(i+b);m++;//System.out.print(i + " ");
			    }
			    i++;
		}
		frequents[m]="o";
		return frequents;
	}
	public char Complement(char a){ //complement symbols
		char c = ' ';
		if(a=='A')c='T';
		if(a=='T')c='A';
		if(a=='G')c='C';
		if(a=='C')c='G';
        return c;
	}
	public String ComplementReverse(String dna){ //complement reversal of string
		String s="";
		int i=dna.length()-1;
	    while(i>=0){
	    	s=s+Complement(dna.charAt(i));
	    	i--;
	    }
		return s;
	}
	public void PatternMatch(String pattern,String dna){ //returns starting indices of pattern occurrences
		int i=0;
		for(i=0;i<(dna.length()-pattern.length());i++){
			if(pattern.equals(dna.substring(i,i+pattern.length())))System.out.print(i +  " ");
		}
	}
	public void FindClump(String genome,int k,int l , int t){ //finds k-mers forming clump
		int i=0,j=0;
		String[] kmers =new String[1000000];
		while(i<genome.length()-l){
			String window = genome.substring(i, i+l);
			String[] temp = new String[10000];
			temp = FrequentPatterns(window, k, t);
			int p=0;
			while(temp[p]!="o"){
				int z=0,h=0;
				while(z<j){
				if(kmers[z].equals(temp[p]))h=1;
				z++;}
				if(h==0){kmers[j]=temp[p];
				j++;}p++;
			}
			i++;
		}
		kmers[j]="o";i=0;
		while(kmers[i]!="o"){
			if(kmers[i]!=null)System.out.print(kmers[i]+" ");
			i++;
		}
	}
    public int[] Skew(String dna){    //#G - #C array or skew array
    	int[] skew = new int[dna.length()+1];
    	int i=0;skew[i]=0;i=1;
 
    	while(i<=dna.length()){
    		if(dna.charAt(i-1)=='G')skew[i]=skew[i-1]+1;
    		else if(dna.charAt(i-1)=='C')skew[i]=skew[i-1]-1;
    		else skew[i]=skew[i-1];
    		i++;
    	}
    	return skew;
    }
    public int HammingDistance(String a,String b){ // distance among strings
    	int ham=0;
    	for(int i =0;i<a.length();i++){
    		if(a.charAt(i)!=b.charAt(i))ham++;
    	}
    	return ham;
    }
    public int[] SkewMin(int[] skew){ //minimum skew value positions of dna
    	int[] min = new int[100];
    	min[0]=0;int tempmin = skew[0];int count =1;
    	for(int i=1;i<skew.length;i++){
    		int temp = skew[i];
    		if(temp==tempmin){
    			min[count]= i;count++;
    		}
    		else if(temp<tempmin){
    			tempmin = temp;
    			for(int j=0;j<count;j++){
    				min[j]=0;
    			}
    			count=1;
    			min[0]=i;
    		}
    	}
    		min[count]=-1;
    	return min;
    }
    public int[] ApproxMatchIndices(String pattern,String dna,int d){//starting indices approx matches of k-mer pattern
    	int[] ind = new int[10000];int count = 0;
    	for(int i=0;i<dna.length()-pattern.length()+1;i++){
    		if(HammingDistance(pattern,dna.substring(i,i+pattern.length()))<=d){
    			ind[count]=i;
    			count++;
    		}
    	}ind[count]=-1;
    	return ind;
    }
    public int ApproxMatchCount(String pattern,String dna,int d){
    	int[] arr =  ApproxMatchIndices(pattern, dna, d);
    	int c = 0;
    	while(arr[c]!=-1)c++;
    	return c;
    }
    
}
