package ac.jp.titech.ntt;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

import org.knowceans.corpus.LabelNumCorpus;
import org.knowceans.corpus.NumCorpus;
import org.knowceans.topics.cgen.Pam5GibbsSampler;
import org.knowceans.util.CokusRandom;
import org.knowceans.util.StopWatch;

public class RunPam5GibbsSampleCrossValidation {
	public static void main(String args[]) throws IOException {
		SProperties props = new SProperties();
		int fold, K, L, J;
		double perplexity[];
		String corpusFile, perplexityFile;

		props.load(args);

		props.checkNullKeys(new String[] { "copus_file", "fold", "K", "L", "J", "perplexity_file" });
		corpusFile = props.getProperty("copus_file");
		fold = Integer.valueOf(props.getProperty("fold"));
		K = Integer.valueOf(props.getProperty("K"));
		L = Integer.valueOf(props.getProperty("L"));
		J = Integer.valueOf(props.getProperty("J"));
		perplexityFile = props.getProperty("perplexity_file");
		
		perplexity = new double[fold + 1];

		// parameters
		double alpha = 0.01;
		double gamma = 0.01;
		double delta = 0.01;
		double beta = 0.01;

		int niter = 2000, niterq = 2000, V, w[][], wq[][];

		// set up corpus
		Random rand = new CokusRandom();
		NumCorpus corpus, train, test;

		corpus = new NumCorpus(corpusFile);

		for (int ii = 0; ii <= fold; ii++) {
			if (ii == 0) {
				corpus.split(fold, 0, rand);
			} else {
				corpus.split(fold, ii, null);
			}

			train = (NumCorpus) corpus.getTrainCorpus();
			test = (NumCorpus) corpus.getTestCorpus();

			// TODO: adjust data source
			w = train.getDocWords(rand);
			wq = test.getDocWords(rand);
			V = corpus.getNumTerms();

			// create sampler
			Pam5GibbsSampler gs = new Pam5GibbsSampler(alpha, gamma, delta, K, L, J, beta, w, wq, V, rand);
			gs.init();
			System.out.println(gs);
			// initial test
			gs.initq();
			gs.runq(niterq);
			System.out.println(gs.ppx());
			// run sampler
			StopWatch.start();
			gs.init();
			gs.run(niter);
			System.out.println(StopWatch.format(StopWatch.stop()));
			// final test
			gs.initq();
			gs.runq(niterq);
			perplexity[ii] = gs.ppx();
			System.out.println(gs.ppx());
			System.out.println(gs);
		}
		
		FileWriter out = new FileWriter(perplexityFile);
		for(int ii=0; ii<fold; ii++) {
			out.write(perplexity[ii]+" ");
		}
		out.close();
	}
}
