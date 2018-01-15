package machinelearning;

public class testfor {

	// Computes log Prob(observation|model) using the Forward algorithm
	// Problem 1 Given the observation sequence O and model, calculate the P(O|model)
	// Forwardscale method
	public static void main(String[] args) throws Exception {

		if(args.length!=2)
		{
			System.out.println("Usage error");
			System.out.println("Usage: java -jar I-HMM_testvit.jar <model.hmm> <obs.seq>\n");
			System.exit(0);
		}
		
//		HMM hmm = new HMM("D:/data/t2.hmm");	
		HMM hmm = new HMM(args[0]);
//		hmm.Output_HMMMode();		
		//hmm.Output_HMMMode();		
//		int[] obs = ReadSequence("E:/t2.1500.seq");
		int[] obs = MHMM.ReadSequence(args[1]);
		
//		System.out.println(obs.length);
		
		int[] q = new int[obs.length];
		double[][] delta = new double[obs.length+1][hmm.getN()+1];
		int[][] psi = new int[obs.length+1][hmm.getN()+1];			
		double[][] alpha = new double[obs.length][hmm.getN()+1];
		double[] scale = new double[obs.length];				
		// Question 3		
		// gamma = dmatrix(1, T, 1, hmm.N);
		double[][] gamma = new double[obs.length][hmm.getN()+1];
		double[][] beta = new double[obs.length][hmm.getN()+1];		
		String[] pred = MHMM.BaumWelch(hmm,obs.length-1,obs,alpha, beta,gamma).split("\\s+");		
//		System.out.println("Training model:");
//		hmm.Output_HMMModel();		
		///// final log score divid the length of the sequence would be better
		System.out.println("\nIteration=" + pred[0] + "initial Log" + pred[1] + " final Log" + pred[2]);
		
	}
}
