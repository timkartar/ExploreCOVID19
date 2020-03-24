import java.util.Random;

public class MotifSearch {
	Random random = new Random();
	static OriSearch ori = new OriSearch();
	
	public String[] Neighbours(String Pattern,int d){
		String[] neighbours = new String[1000000];
		if(d == 0){
			neighbours[0]=Pattern;neighbours[1]="o";return neighbours;
		}
		if(Pattern.length()==1){
			neighbours[0]="A";neighbours[1]="C";
			neighbours[2]="G";neighbours[3]="T";
			neighbours[4]="o";
			return neighbours;
		}
		String[] suffix_neighbours = Neighbours(Pattern.substring(1),d);
		int i=0,count=0;
		while(suffix_neighbours[i]!="o"){
			if(ori.HammingDistance(suffix_neighbours[i],Pattern.substring(1))<d){
				for(int j=0;j<4;j++){
					neighbours[count]=ori.DigtoChar(j) + suffix_neighbours[i];count++;
				}
			}
			else{
				neighbours[count]=Pattern.charAt(0) + suffix_neighbours[i];count++;
			}
			i++;
		}
		neighbours[count]="o";
		return neighbours;
	}
    public String[] MotifEnumeration(String[] dna,int k,int d){
    	/*System.out.println("checking");*/
    	 String[] Patterns = new String[10000];
    	 int i=0,count=0;
    	 int l=0;int m=0;
    	 while(l<dna.length){/*System.out.println("checking1");*/
    	   while(i<dna[l].length()-k+1){/*System.out.println("checking2");*/
    		 String[] d_match=Neighbours(dna[l].substring(i,i+k),d);
    		 int j=0;
    		 while(d_match[j]!="o"){/*System.out.println("checking3");*/
    			   for(int g=0;g<dna.length;g++){/*System.out.println("checking4");*/
    				   int f=ori.ApproxMatchCount(d_match[j],dna[g], d);//System.out.print(d_match[j] + " ");
    			       if(f==0){//System.out.print(f + " ");
    				      m=1;break;
    			       }
    		        }//System.out.print(m);
    		       if(m==0){Patterns[count]=d_match[j];count++;//System.out.println("checking5 "+count);
    			       }m=0;
    		         j++;//System.out.println("checking6");
    			    }
    		 i++;//System.out.println("checking7");
    	    }
    	   l++;//System.out.println("checking8");
    	  }
    	 Patterns[count]="o";
    	 return Patterns;
    }
    public String[] RemoveRepeats(String[] dna){
    	String[] str2 = new String[1000];
    	int i=0,count=0,m=0;str2[0]=dna[0];i++;count++;
    	while(dna[i]!="o"){
    		for(int j=0;j<count;j++){
    			if(str2[j].equals(dna[i]))m=1;
    		}
    		if(m==0){
    			str2[count]=dna[i];count++;
    		}
    	    i++;
    	}str2[count]="o";
    	return str2;
    }
    public int DisPatStrings(String pat,String[] dna){
    	int dis = 0,i=0;
    	int k = pat.length();
    	while(i<dna.length){
    		int ham=k;
    		int j=0;
    		for(j=0;j<dna[i].length()-k+1;j++){
    			int h = ori.HammingDistance(pat, dna[i].substring(j,j+k));
    			if(ham>h){
    				ham=h;
    			}
    		 }
    		dis = dis + ham;i++;
    	}
    	return dis;
    }
    public String MedianString(String[] dna,int k){
    	String pat =" ",med=" ";
    	int distance = k;double a=Math.pow(4, k);
    	for(double i=a;i<a+a;i++){
    		pat=ori.NumbertoString(i);
    		//System.out.print(pat+" ");
    		int dis=DisPatStrings(pat,dna);
    		if(distance>dis){distance=dis;med=pat;}
    	}
    	return med;
    }
    public char DigtoChar(int a){
    	char c='A';
    	if(a==1)c='C';
    	if(a==2)c='G';
    	if(a==3)c='T';
    	return c;
    }
    public int ChartoDig(char a){
    	int i=0;
    	if(a=='A')i=0;
    	if(a=='C')i=1;
    	if(a=='G')i=2;
    	if(a=='T')i=3;
    	if(a=='o')i=4;
    	//System.out.print(i);
    	return i;
    }
    public double ProfileProbability(String dna,double[][] profile){
    	double prb =1;
    	for(int i=0;i<dna.length();i++){
    		prb=prb*profile[ChartoDig(dna.charAt(i))][i];//System.out.print(profile[ChartoDig(dna.charAt(i))][i]+" ");
    	}
    	return prb;
    }
    public String ProfileMostKmer(String dna,int k,double[][] profile){
       	String kmer=dna.substring(0,k);
    	int i=0;
    	double prob=0;
    	for(i=0;i<dna.length()-k+1;i++){
    		double tempProb=ProfileProbability(dna.substring(i,i+k),profile);
    		//System.out.print(tempProb);//System.out.print(profile[0][1]);
    		if(tempProb>prob){kmer = dna.substring(i, i+k);prob=tempProb;}
    	}
    	
    	return kmer;
    }
	public int Score(String[] motifs,int k){
		int score=0;
		String con="" ;
		double[][] profile = MakeProfile(motifs,k);
		/*for(int i=0;i<4;i++){
			for(int j=0;j<k;j++){
				System.out.print(profile[i][j]+" ");
			}System.out.print("\n");
		}*/
		for(int i=0;i<k;i++){
			double temp=0;
			int dig=0;
			for(int j=0;j<4;j++){
				if(profile[j][i]>temp){
					temp=profile[j][i];//System.out.print(temp);
					dig=j;
				}
			}
			con=con+DigtoChar(dig);
		}
		//System.out.print(con+" ");
		int i=0;
		while(motifs[i]!="o"){
			i++;
		}i=i-1;
		for(int j=0;j<i;j++){
			//System.out.print(j+" ");
			score=score+ori.HammingDistance(con,motifs[j]);
		}
		return score;
	}
	public double[][] MakeProfile(String[] motifs,int k){
		double[][] profile = new double[4][k];
		int l=0;
		while(motifs[l]!="o"){
			l++;
		}//System.out.print(l+" ");
		for(int i=0;i<k;i++){
			for(int j=0;j<l;j++){
				//System.out.print(/*ChartoDig(motifs[j].charAt(i))*/i+" ");
				
				profile[ChartoDig((motifs[j]).charAt(i))][i]=profile[ChartoDig(motifs[j].charAt(i))][i]+1;
				
			}//System.out.print(" ");
		}
		for(int i=0;i<4;i++){
			for(int j=0;j<k;j++){
				profile[i][j]++;
				profile[i][j]=profile[i][j]/(l+4);
			}
		}
		return profile;
	}
    public String[] GreedyMotifSearch(String[]dna,int k,int t){
		String[] BestMotifs = new String[t+1];
		for(int i=0;i<t;i++){
		 BestMotifs[i]=dna[i].substring(0,k);
		}BestMotifs[t]="o";
	    String[] Motifs = new String[t+1];
	    for(int i=0;i<(dna[0].length()-k+1);i++){
	    	Motifs[0]=dna[0].substring(i,i+k);
	    	Motifs[1]="o";
	    	
	    	for(int j=1;j<t;j++){
	    		double[][] profile = MakeProfile(Motifs,k);
	    		for(int l=0;l<4;l++){
	    		for(int m=0;m<k;m++){
	    			//System.out.print(profile[l][m]+" ");
	    		}//System.out.print("\n");
	    	 }
	    		//System.out.print(ProfileMostKmer(dna[j],k,profile)+"avada");
	    		Motifs[j]=ProfileMostKmer(dna[j],k,profile);
	    		//System.out.print(" "+ Motifs[j]);
	    		Motifs[j+1]="o";
	    	} // double[][] profile=MakeProfile(Motifs,k);
	    	//System.out.print(Score(Motifs)+" ");
	    	//System.out.print(Score(BestMotifs)+" ");
	    	if(Score(Motifs,k)<Score(BestMotifs,k)){
	    		//System.out.print("yes");
	    		int p=0;
	    		while(Motifs[p]!="o"){
	    			BestMotifs[p]=Motifs[p];
	    			p++;
	    		}
	    	}
	    }
	    return BestMotifs;
	}
    public String[] Motifs(double[][] profile,String[] dna,int k){
    	String[] motifs = new String[dna.length];
    	for(int i=0;i<dna.length-1;i++){
    		motifs[i]=ProfileMostKmer(dna[i],k,profile);
    	}
    	motifs[dna.length-1]="o";
    	return motifs;
    	}
    public String[] RandomMotifSearch(String[] dna,int k,int t){
    String[] BestMotifs =  new String[t+1];
    String[] motifs = new String[t+1];
    for(int i=0;i<t;i++){
    	int r = random.nextInt(dna[0].length()-k+1);
    	motifs[i]=dna[i].substring(r,r+k);
    }
    motifs[t]="o";
    int p=0;
	while(p<t+1){
		BestMotifs[p]=motifs[p];
		p++;
	}
	double[][] profile = new double[4][k];
	while(true){
		profile=MakeProfile(motifs,k);
		motifs=Motifs(profile,dna,k);
		if(Score(motifs,k)<Score(BestMotifs,k)){
			p=0;
			while(p<t+1){
				BestMotifs[p]=motifs[p];
				p++;
			}
		}
		else
			return BestMotifs;
	}
    
   }
}