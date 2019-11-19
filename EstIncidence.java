package mnbprojects;
import java.util.*;
import java.io.*;


/**
 * Utility for estimating disease incidence given a set of 'true' minor allele frequencies of causative alleles.
 * This program takes, from the command line, a random seed, a file with minor allele frequencies 
 * (one per line, in any acceptable double format), and the total number of people wanted in your simulated 
 * cohort (50-100 typically).
 * 
 * Requires: cohen_biometrics_1960.txt to be in same directory.
 * 
 * Given a set of pathogenic allele frequencies, this program will 
 * 1) generate a "gnomad-like" cohort of unaffected individuals
 * 2) Generate a cohort of affected individuals
 * 3) use this data to estimate the total pathogenic allele frequency 
 * 
 * 
 * @author mbainbridge
 *
 */
public class EstIncidence 
{
	public static double mafcutoffforrare = 1e-5; //this is used for the lambda calculation.  should be set to something close to the limits of detection for gnomad
	public static Random rng; // = new Random(123454321);
	public static void main(String[] args) throws Exception
	{
		if(args.length != 3)
		{
			System.err.println("Usage: <RndSeed> <FileOfMAFS> <total people wanted>");
			System.exit(1);
		}
		rng = new Random(Integer.parseInt(args[0]));
		
		int totalPeopleWanted = Integer.parseInt(args[2]);
		int totalAffected = 0;
		int gnoamdsize = 600000;
		int tries = 0;
		double[] realInc ;
		realInc = genRealInc();
		double sum = 0;
		double[] cumudistn = new double[realInc.length];
		for(int i = 0; i < realInc.length; i++)
		{
			sum+=realInc[i];
			cumudistn[i] = sum;
		}
		
		//generate the gnomad-like database
		double[] gnomad = EstIncidence.genGnomad(realInc,  gnoamdsize, cumudistn);
		
		//generate the cohort
		int[] simHitCnt = new int[realInc.length];
		while(totalAffected < totalPeopleWanted)
		{
			tries++;
			if(tries%1000000 == 0) System.out.print(".");
			int index1 = hits(cumudistn);
			int index2 = hits(cumudistn);
			if(index1 != -1 && index2 != -1)
			{
				System.out.print("!");
				totalAffected++;
				simHitCnt[index1]++;
				simHitCnt[index2]++;	
			}
		}
		//print out, for each allele, the real MAF, the MAF in our simulated gnomad, and the numebr of individuals in oru affected cohort with that allele.
		for(int i = 0; i < realInc.length; i++)
		{
			System.out.println(i+" "+realInc[i]+" "+gnomad[i] +" "+simHitCnt[i]);
		}
		System.out.println("\nTotal tries: "+tries); //the number of people samples until we had our cohort
		
		double totMAF = 0;
		int zeroMAFcnts  = 0;
		int MAFcnts = 0;
		int uniqueobsv = 0;
		int zeroobs = 0; 
		double realObservedMAF = 0.0d;  //real or actual observed MAF as opoosed to our estimator
		double realUnObservedMAF = 0.0d;
		for(int i = 0; i < realInc.length; i++)
		{
			if(simHitCnt[i] != 0)
			{
				uniqueobsv++;
				realObservedMAF+=realInc[i];
				if(gnomad[i] == 0)
				{
					zeroMAFcnts+=simHitCnt[i];
				}
				else
				{
					totMAF+=gnomad[i];
					MAFcnts+=simHitCnt[i];
				}		
			}	
			else
			{
				realUnObservedMAF+=realInc[i];
				zeroobs++;
			}
		}
		/*
		 * Calculate the missing maf for observed alleles by making it proportional to the observed alleles with known maf
		 */
		double missingMAF = totMAF*zeroMAFcnts/MAFcnts;
		double adjMAF = missingMAF+totMAF;
		double observedDiff = realObservedMAF - adjMAF;
		System.out.println("Total real MAF: "+sum);
		System.out.println("Actual_observedMAF: "+realObservedMAF);
		System.out.println("Estimated MAF: "+adjMAF);
		System.out.println("Diff_between_estimatedObs_and_real: "+observedDiff);
		double pcdiff = observedDiff/adjMAF;
		System.out.println("PC_Diff_between_estimatedObs_and_real: "+pcdiff);
		//double LL = lambdaEstimator(gnomad, simHitCnt);
		double L = lambdaEstimator2(simHitCnt, gnomad);
		//System.out.println("Lambda: "+L + " Old estimator "+ LL);
		double k0 = Math.pow(Math.E,L*-1.0d);
		double propzero = (uniqueobsv)/(1.0d-k0)*k0;
		System.out.println("Number of missing alleles: "+propzero + "  actualy: "+zeroobs);
		double diffMissingAlleles = zeroobs - propzero;
		pcdiff = diffMissingAlleles/(double)zeroobs;
		System.out.println("Diff_missing_alleles: "+diffMissingAlleles);
		System.out.println("PC_Diff_missing_alleles: "+pcdiff);
		double MAFforMissingValues = mafformissingvalues(simHitCnt,gnomad);
		System.out.println("MAF for missing values: "+MAFforMissingValues);
		double additionalMAF = propzero*MAFforMissingValues;
		System.out.println("Add MAF for missing: "+additionalMAF);
		System.out.println("Actually_unobservedMAF: "+realUnObservedMAF);
		double diffMissingMAF = realUnObservedMAF - additionalMAF;
		System.out.println("Diff_maf_missing: "+diffMissingMAF);
		pcdiff = diffMissingMAF/realUnObservedMAF;
		System.out.println("PC_Diff_maf_missing: "+pcdiff);
		double finMAF = adjMAF+additionalMAF;
		System.out.println("Total final MAF: "+finMAF);
	}
	
	static double[][] lambdaME = new double[0][0]; 
	
	/**
	 * Estimates the most likely lamdba from a truncated poisson distribution.  Cohen, Biometrics, June 1960.
	 * @param cnts  distribution of allele observations, i.e how many alleles were observed once, twice, three times, etc.
	 * @param gnomad  Our generated gnomad cohort
	 * @return Lambda parameter for Poisson distribution
	 * @throws Exception
	 */
	public static double lambdaEstimator2(int[] cnts, double[] gnomad) throws Exception
	{
		//Load the precomputed lambdas from cohen 1960 if not already laoded.
		if(lambdaME.length == 0)
		{
			BufferedReader br = new BufferedReader(new FileReader("cohen_biometrics_1960.txt"));
			String line;
			int cnt = 0;
			while((line = br.readLine())!=null)
			{
				cnt++;
			}
			br.close();
			br = new BufferedReader(new FileReader("cohen_biometrics_1960.txt"));
			lambdaME = new double[cnt][2];
			cnt = 0;
			while((line = br.readLine())!=null)
			{
				String[] s = line.split("[ \t]+");
				lambdaME[cnt][0] = Double.parseDouble(s[0]);
				lambdaME[cnt][1] = Double.parseDouble(s[1]);
				cnt++;
			}
		}
		int sum = 0;
		int non0cnt = 0;
		//generate teh average number of counts but only for alleles below our 'rare' threshold.
		for(int i = 0; i < cnts.length; i++)
		{
			if(gnomad[i] < mafcutoffforrare)
			{
				sum+=cnts[i];	
				if(cnts[i] > 0) non0cnt++;
			}
		}
		double avg = (double)sum/(double)non0cnt;
		double min = 100000.0d;
		
		System.out.println("lambda2.avg: "+avg);
		
		//find the closest lambda from the table
		double bestlambda = -1.0d;
		for(int i = 0; i < lambdaME.length; i++)
		{
			double xxx = Math.abs(avg-lambdaME[i][0]);
			if(xxx < min)
			{
				min = xxx;
				bestlambda = lambdaME[i][1];
				//System.out.println("bestlambda: "+bestlambda+" i:"+i+" "+lambdaME[i][0]+" "+lambdaME[i][1]);
			}
			
		}
		return bestlambda;	
	}
	
	/**
	 * Calculates the MAF value that should be used for missing alleles
	 * @param cnts  Simulated hit count for each allele
	 * @param gnomad  Simualetd gnoamd-liek cohort
	 * @param unobsalleles 
	 * @return
	 */
	public static double mafformissingvalues(int[] cnts, double[] gnomad)
	{
		int sumcnts = 0;
		double summafs = 0;
		for(int i = 0; i < cnts.length; i++)
		{
			//if(cnts[i] < 2)
			if(gnomad[i]<mafcutoffforrare)
			{
				summafs+=cnts[i]*gnomad[i];
				sumcnts+=cnts[i];
			}
		}
		double ds = (double) sumcnts;		
		System.out.println("Sum mafs and sum cnt "+summafs+" "+ds);
		return summafs/ds; 
	}
	
	
	/**
	 * Loads MAFs from a text file		
	 * @param fname File with MAF data
	 * @return A set of MAFs
	 * @throws Exception
	 */
	public static double[] loadRfromFile(String fname) throws Exception
	{
		BufferedReader br = new BufferedReader(new FileReader(fname));
		String line;
		int cnt = 0;
		while((line = br.readLine())!= null)
		{
			cnt++;			
		}
		br.close();
		double[] D = new double[cnt];
		br = new BufferedReader(new FileReader(fname));
		cnt = 0;
		while((line = br.readLine())!= null)
		{
			double d = Double.parseDouble(line);
			D[cnt]=d;
			cnt++;
		}
		return D;
	}
	
	 /**
	  * Determines whether a person has a pathogenic allele or not.
	  * 
	  * @param cumudistn The cumulative distribution of pathogenic MAFs
	  * @return The allele number where the hit occurred or -1 if no hits
	  */
	public static int hits( double[] cumudistn)
	{
		double rnd = rng.nextDouble();
		for(int i = 0; i < cumudistn.length; i++)
		{
			
			if(rnd <= cumudistn[i]) return i ;
		}
		return -1;	
	}
	
	
	/**
	 * @deprecated
	 * @param gnomad
	 * @param cnts
	 * @return
	 */
	public static double lambdaEstimator(double[] gnomad, int[] cnts)
	{
		int twopluscnt = 0;
		int onecnt = 0;
		boolean debug = true;
		for(int i = 0; i < gnomad.length; i++)
		{
			if(gnomad[i] < mafcutoffforrare)
			{
				if(cnts[i] > 1) twopluscnt++;
				if(cnts[i] == 1) onecnt++;
			}
		}
		
		double ratio = (double)onecnt/(double)twopluscnt;
		if(debug)
		{
			System.out.println("one cnt and 2 cnt rat :"+onecnt+"  "+twopluscnt+ " "+ratio);
		}
		double min = 10000;
		double bestlambda = 0;
		for(int i = 1; i < 200; i++ )
		{
			double lambda = (double) i/100.0d;
			double k1 = Math.pow(Math.E,(-1.0d*lambda))*Math.pow(lambda,1.0d)/(double)EstIncidence.factorial(1);
			double sumk2345 = 0;
			for(int j = 2; j < 6; j++)
			{
				sumk2345+=Math.pow(Math.E,(-1.0d*lambda))*Math.pow(lambda,(double)j)/(double)EstIncidence.factorial(j);
			}
			double calcrat = k1/sumk2345;
			if(debug)
			{
				if(i >49 && i < 86)
				{
					System.out.println(i+" "+k1+" "+calcrat);
				}
			}
			double diff = Math.abs(calcrat - ratio);
			//if(diff > min) return (double) (i-1)/100.0d;
			if(diff < min)
			{
				min = diff;
				bestlambda = lambda;
			}
		}
		return bestlambda;
	}
	
	public static long factorial(int number) 
	{
		long result = 1;
        for (int factor = 2; factor <= number; factor++) {
            result *= factor;
        }
        return result;
    }
	
	/**
	 * Generates a simulated gnomad like database from the given MAFs
	 * @param real  Real MAF distribution
	 * @param numberofpeople Number of people wanted in gnomad cohort
	 * @param cumu  The cumulative distribution of pathogenic MAFs
	 * @return
	 * @throws Exception
	 */
	public static double[] genGnomad(double[] real, int numberofpeople, double[] cumu) throws Exception
	{
		int[] hitcnt = new int[real.length];
		double[] gnomad = new double[real.length];
		for(int i = 0; i < numberofpeople; i++)
		{
			for(int j = 0; j < 2; j++)
			{
				int hit = hits(cumu);
				if(hit != -1)
				{
					//System.out.println("1111");
					hitcnt[hit]++;
					break;
				}
			}
		}
		for(int i = 0; i < real.length; i++)
		{
			gnomad[i] = (double)hitcnt[i]/(numberofpeople*2);
		}
		return gnomad;
	}
	
	/**
	 * 
	 * Generates a distribution of real pathogenic MAFs
	 * @deprecated
	 * @return
	 * @throws Exception
	 */
	public static double[] genRealInc() throws Exception
	{
		//generate real incidence
		double[] realInc = new double[200];
		realInc[0] =1.00E-04;
		realInc[1] =5E-05;
		realInc[2] =5E-05;
		realInc[3] =2E-05;
		double bonus = 0.0;
		int cnt = 4;
		while(bonus < 3.9E-4 && cnt < realInc.length)
		{
			double g = rng.nextGaussian();
			
			double value = 3E-6+1.0E-6*g;
			if(value < 0) continue;
			bonus+=value;
			realInc[cnt]=value;
			cnt++;
			
		}
		System.out.println("count_and_bonus: "+cnt+"  "+bonus);
		return realInc;
		
		
	}
}
