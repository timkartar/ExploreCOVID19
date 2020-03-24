import java.util.Scanner;
public class Main {
    static MotifSearch motif = new MotifSearch();
	public static void main(String[] args){
	//System.out.print(motif.a);
	Scanner scanner = new Scanner(System.in);
	int k =scanner.nextInt();
	int t=scanner.nextInt();
	//double[][] profile=new double[4][t];//System.out.print(k);
	String[] dna=new String[t+1];
	for(int i=0;i<t;i++){
	   dna[i]=scanner.next();
	//System.out.print("\n");
}   //System.out.print(profile[0][1]);
	//scanner.close();
	dna[t]="o";
	//double[][] profile=motif.MakeProfile(dna, 12);
	/*for(int i=0;i<4;i++){
		for(int j=0;j<t;j++){
			profile[i][j]=scanner.nextDouble();
		}System.out.print("\n");
	}*/
	scanner.close();
	//System.out.print(dna[0].length());
	//System.out.print(motif.ProfileProbability(dna, profile));
	//String[] str = motif.GreedyMotifSearch(dna, k, t);
	//for(int j=0;str[j]!="o";j++)System.out.println(str[j]);
	String[] BestMotifs = motif.RandomMotifSearch(dna, k, t);
	for(int i=0;i<999;i++){
		String[]Motifs = motif.RandomMotifSearch(dna, k, t);
		if(motif.Score(Motifs,k)<motif.Score(BestMotifs,k)){
			int p=0;//System.out.print("s"+motif.Score(Motifs,k));
			while(p<t){
				
				BestMotifs[p]=Motifs[p];
				p++;
			}
		}
	 }
	for(int j=0;j<t;j++){
		System.out.println(BestMotifs[j]);
	}
	System.out.print(motif.Score(BestMotifs,k));
    } 
}

