/*
 * Pam5.java -- PAM5 Gibbs sampling implementation
 * Copyright (C) 2009-2011 Gregor Heinrich, gregor :: arbylon : net
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License version 3 as
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.knowceans.topics.cgen;

// imports
import java.util.Random;
import org.knowceans.corpus.LabelNumCorpus;
import org.knowceans.util.CokusRandom;
import org.knowceans.util.DirichletEstimation;
import org.knowceans.util.StopWatch;

// 
/**
 * Generated Gibbs sampler for the PAM5 model. 5-level pachinko allocation
 * model, which extends the basic 4-level model by an additional hierarchy
 * level. Testing C1B structures and hyperparameter grouping.
 * 
 * Mixture network specification:
 * 
 * <pre>
 * m >> (theta[m] | alpha) >> x[m][n]
 * m, x[m][n] >> (thetax[m,x] | gamma[x]) >> y[m][n]
 * m, y[m][n] >> (thetay[m,y] | delta[y]) >> z[m][n]
 * z[m][n] >> (phi[z] | beta) >> w[m][n]
 * 
 * 
 * elements:
 * 
 *    document m
 *    type: {ROOT|E3COUPLED}
 *    parents: (root)
 *    children: theta thetax thetay
 *    range: M
 * ----
 *    doc-topic (theta[m] | alpha)
 *    type: {SEQUENCE|C1ROOT}
 *    parents: m
 *    children: x
 *    components: theta, domain: M, range: K
 *    counts: nmx, sum: null
 *    index: m, selector: null
 *    hyperparams: alpha, dimension 1, fixed: false 
 * ----
 *    topic-topic (thetax[m,x] | gamma[x])
 *    type: {SEQUENCE|C1BSEQADD}
 *    parents: m x
 *    children: y
 *    components: thetax, domain: M *  K, range: L
 *    counts: nmxy, sum: nmxysum
 *    index: mx, selector: m,x
 *    hyperparams: gamma, dimension: K * L, selector: x, fixed: false 
 * ----
 *    topic-topic (thetay[m,y] | delta[y])
 *    type: {SEQUENCE|C1BSEQADD}
 *    parents: m y
 *    children: z
 *    components: thetay, domain: M *  L, range: J
 *    counts: nmyz, sum: nmyzsum
 *    index: my, selector: m,y
 *    hyperparams: delta, dimension: L * J, selector: y, fixed: false 
 * ----
 *    topic x[m][n]
 *    type: {HIDDEN|E1SINGLE}
 *    parents: theta
 *    children: thetax
 *    range: K
 * ----
 *    topic y[m][n]
 *    type: {HIDDEN|E1SINGLE}
 *    parents: thetax
 *    children: thetay
 *    range: L
 * ----
 *    topic z[m][n]
 *    type: {HIDDEN|E1SINGLE}
 *    parents: thetay
 *    children: phi
 *    range: J
 * ----
 *    topic-word (phi[z] | beta)
 *    type: {TOPIC|C1ASINGLE}
 *    parents: z
 *    children: w
 *    components: phi, domain: J, range: V
 *    counts: nzw, sum: nzwsum
 *    index: z, selector: null
 *    hyperparams: beta, dimension 1, fixed: false 
 * ----
 *    word w[m][n]
 *    type: {VISIBLE|E1SINGLE}
 *    parents: phi
 *    children: (leaf)
 *    range: V
 * ----
 * sequences:
 * 
 * words [m][n]
 * parent: (root), children: []
 * edges: m x y z w
 * </pre>
 * 
 * @author Gregor Heinrich, gregor :: arbylon : net (via MixNetKernelGenerator)
 * @date generated on 11 Mar 2011
 */
// constructor
public class Pam5GibbsSampler {

	// //////////////// fields //////////////////
	// fields
	Random rand;
	int iter;
	int niter;
	// root edge: document
	int M;
	int Mq;
	// sequence node: doc-topic
	int[][] nmx;
	int[][] nmxq;
	double alpha;
	double alphasum;// sequence node: topic-topic
	int[][] nmxy;
	int[][] nmxyq;
	int[] nmxysum;
	int[] nmxysumq;
	double[][] gamma;
	double[] gammasum;// sequence node: topic-topic
	int[][] nmyz;
	int[][] nmyzq;
	int[] nmyzsum;
	int[] nmyzsumq;
	double[][] delta;
	double[] deltasum;// hidden edge: topic
	int[][] x;
	int[][] xq;
	int K;
	// hidden edge: topic
	int[][] y;
	int[][] yq;
	int L;
	// hidden edge: topic
	int[][] z;
	int[][] zq;
	int J;
	// topic node: topic-word
	int[][] nzw;
	int[] nzwsum;
	double beta;
	double betasum;
	double[][] phi;
	// visible edge: word
	int[][] w;
	int[][] wq;
	int V;
	// sequence: words
	int W;
	int Wq;
	// sampling weights
	double[][][] pp;

	// //////////////// main //////////////////
	// standard main routine
	public static void main(String[] args) {
		String filebase = "nips/nips";
		Random rand = new CokusRandom();
		// set up corpus
		LabelNumCorpus corpus = new LabelNumCorpus(filebase);
		corpus.split(10, 1, rand);
		LabelNumCorpus train = (LabelNumCorpus) corpus.getTrainCorpus();
		LabelNumCorpus test = (LabelNumCorpus) corpus.getTestCorpus();
		// TODO: adjust data source
		int[][] w = train.getDocWords(rand);
		int[][] wq = test.getDocWords(rand);
		int V = corpus.getNumTerms();
		// parameters
		int K = 20;
		int L = 20;
		int J = 20;
		double alpha = 0.1;
		double gamma = 0.1;
		double delta = 0.1;
		double beta = 0.1;

		int niter = 10, niterq = 5;
		// create sampler
		Pam5GibbsSampler gs = new Pam5GibbsSampler(alpha, gamma, delta, K, L,
				J, beta, w, wq, V, rand);
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
		System.out.println(gs.ppx());
		System.out.println(gs);

	} // main

	// //////////////// constructor //////////////////
	public Pam5GibbsSampler(double alpha, double gamma, double delta, int K,
			int L, int J, double beta, int[][] w, int[][] wq, int V, Random rand) {
		// assign
		this.K = K;
		this.L = L;
		this.J = J;
		this.w = w;
		this.wq = wq;
		this.V = V;
		this.alpha = alpha;
		this.alphasum = K * alpha;
		this.gamma = new double[K][L];
		this.gammasum = new double[K];
		for (int mxj = 0; mxj < K; mxj++) {
			for (int t = 0; t < L; t++) {
				this.gamma[mxj][t] = gamma;
			} // for t
			this.gammasum[mxj] = L * gamma;
		} // for mxj
		this.delta = new double[L][J];
		this.deltasum = new double[L];
		for (int myj = 0; myj < L; myj++) {
			for (int t = 0; t < J; t++) {
				this.delta[myj][t] = delta;
			} // for t
			this.deltasum[myj] = J * delta;
		} // for myj
		this.beta = beta;
		this.betasum = V * beta;
		// count tokens
		M = w.length;
		W = 0;
		for (int m = 0; m < M; m++) {
			W += w[m].length;
		}
		Mq = wq.length;
		Wq = 0;
		for (int m = 0; m < Mq; m++) {
			Wq += wq[m].length;
		}
		this.rand = rand;
		// allocate sampling weights
		pp = new double[K][L][J];
	} // c'tor

	// //////////////// initialisation //////////////////
	// initialisation
	public void init() {
		// component selectors
		int mxsel = -1;
		int mxjsel = -1;
		int mysel = -1;
		int myjsel = -1;
		// sequence node
		nmx = new int[M][K];
		// sequence node
		nmxy = new int[M * K][L];
		nmxysum = new int[M * K];
		// sequence node
		nmyz = new int[M * L][J];
		nmyzsum = new int[M * L];
		// hidden edge
		x = new int[M][];
		for (int m = 0; m < M; m++) {
			x[m] = new int[w[m].length];
		}
		// hidden edge
		y = new int[M][];
		for (int m = 0; m < M; m++) {
			y[m] = new int[w[m].length];
		}
		// hidden edge
		z = new int[M][];
		for (int m = 0; m < M; m++) {
			z[m] = new int[w[m].length];
		}
		// topic node
		nzw = new int[J][V];
		nzwsum = new int[J];
		// initialise randomly
		// major loop, sequence [m][n]
		for (int m = 0; m < M; m++) {
			// minor loop, sequence [m][n]
			for (int n = 0; n < w[m].length; n++) {
				// sample edge values
				int hx = rand.nextInt(K);
				int hy = rand.nextInt(L);
				int hz = rand.nextInt(J);
				// assign topics
				x[m][n] = hx;
				y[m][n] = hy;
				z[m][n] = hz;
				// increment counts
				nmx[m][hx]++;
				mxsel = K * m + hx;
				nmxy[mxsel][hy]++;
				nmxysum[mxsel]++;
				mysel = L * m + hy;
				nmyz[mysel][hz]++;
				nmyzsum[mysel]++;
				nzw[hz][w[m][n]]++;
				nzwsum[hz]++;
			} // for n
		} // for m
	} // init

	// //////////////// query initialisation //////////////////
	// initialisation
	public void initq() {
		// component selectors
		int mxsel = -1;
		int mxjsel = -1;
		int mysel = -1;
		int myjsel = -1;
		// sequence node
		nmxq = new int[Mq][K];
		// sequence node
		nmxyq = new int[Mq * K][L];
		nmxysumq = new int[Mq * K];
		// sequence node
		nmyzq = new int[Mq * L][J];
		nmyzsumq = new int[Mq * L];
		// hidden edge
		xq = new int[Mq][];
		for (int m = 0; m < Mq; m++) {
			xq[m] = new int[wq[m].length];
		}
		// hidden edge
		yq = new int[Mq][];
		for (int m = 0; m < Mq; m++) {
			yq[m] = new int[wq[m].length];
		}
		// hidden edge
		zq = new int[Mq][];
		for (int m = 0; m < Mq; m++) {
			zq[m] = new int[wq[m].length];
		}
		// topic node
		// compute parameters
		phi = new double[J][V];
		for (int hz = 0; hz < J; hz++) {
			for (int t = 0; t < V; t++) {
				phi[hz][t] = (nzw[hz][t] + beta) / (nzwsum[hz] + betasum);
			} // t
		} // h
			// initialise randomly
		// major loop, sequence [m][n]
		for (int m = 0; m < Mq; m++) {
			// minor loop, sequence [m][n]
			for (int n = 0; n < wq[m].length; n++) {
				// sample edge values
				int hx = rand.nextInt(K);
				int hy = rand.nextInt(L);
				int hz = rand.nextInt(J);
				// assign topics
				xq[m][n] = hx;
				yq[m][n] = hy;
				zq[m][n] = hz;
				// increment counts
				nmxq[m][hx]++;
				mxsel = K * m + hx;
				nmxyq[mxsel][hy]++;
				nmxysumq[mxsel]++;
				mysel = L * m + hy;
				nmyzq[mysel][hz]++;
				nmyzsumq[mysel]++;

			} // for n
		} // for m
	} // initq

	// //////////////// main kernel //////////////////
	// Gibbs kernel
	public void run(int niter) {
		// iteration loop
		for (int iter = 0; iter < niter; iter++) {
			System.out.println(iter);
			// major loop, sequence [m][n]
			for (int m = 0; m < M; m++) {
				// component selectors
				int mxsel = -1;
				int mxjsel = -1;
				int mysel = -1;
				int myjsel = -1;
				// minor loop, sequence [m][n]
				for (int n = 0; n < w[m].length; n++) {
					double psum;
					double u;
					// assign topics
					int hx = x[m][n];
					int hy = y[m][n];
					int hz = z[m][n];
					// decrement counts
					nmx[m][hx]--;
					mxsel = K * m + hx;
					nmxy[mxsel][hy]--;
					nmxysum[mxsel]--;
					mysel = L * m + hy;
					nmyz[mysel][hz]--;
					nmyzsum[mysel]--;
					nzw[hz][w[m][n]]--;
					nzwsum[hz]--;
					// compute weights
					/*
					 * &p(x_{m,n} \eq x, y_{m,n} \eq y, z_{m,n} \eq z\; |\;\vec
					 * x_{-m,n}, \vec y_{-m,n}, \vec z_{-m,n}, \vec w, \cdot)
					 * \notag\\ &\qquad\propto (n^{-mn}_{m,x} + \alpha ) \cdot
					 * \frac{n^{-mn}_{mx,y} + \gamma_{x, y} }{\sum_{y}
					 * n^{-mn}_{mx,y} + \gamma_{x, y}} \cdot
					 * \frac{n^{-mn}_{my,z} + \delta_{y, z} }{\sum_{z}
					 * n^{-mn}_{my,z} + \delta_{y, z}} \cdot \frac{n^{-mn}_{z,w}
					 * + \beta }{\sum_{w} n^{-mn}_{z,w} + \beta}
					 */
					psum = 0;
					// hidden edge
					for (hx = 0; hx < K; hx++) {
						// hidden edge
						for (hy = 0; hy < L; hy++) {
							// hidden edge
							for (hz = 0; hz < J; hz++) {
								mxsel = K * m + hx;
								mxjsel = hx;
								mysel = L * m + hy;
								myjsel = hy;

								pp[hx][hy][hz] = (nmx[m][hx] + alpha)
										* (nmxy[mxsel][hy] + gamma[mxjsel][hy])
										/ (nmxysum[mxsel] + gammasum[mxjsel])
										* (nmyz[mysel][hz] + delta[myjsel][hz])
										/ (nmyzsum[mysel] + deltasum[myjsel])
										* (nzw[hz][w[m][n]] + beta)
										/ (nzwsum[hz] + betasum);
								psum += pp[hx][hy][hz];
							} // for h
						} // for h
					} // for h

					// sample topics
					u = rand.nextDouble() * psum;
					psum = 0;
					SAMPLED:
					// each edge value
					for (hx = 0; hx < K; hx++) {
						// each edge value
						for (hy = 0; hy < L; hy++) {
							// each edge value
							for (hz = 0; hz < J; hz++) {
								psum += pp[hx][hy][hz];
								if (u <= psum)
									break SAMPLED;
							} // h
						} // h
					} // h

					// assign topics
					x[m][n] = hx;
					y[m][n] = hy;
					z[m][n] = hz;
					// increment counts
					nmx[m][hx]++;
					mxsel = K * m + hx;
					nmxy[mxsel][hy]++;
					nmxysum[mxsel]++;
					mysel = L * m + hy;
					nmyz[mysel][hz]++;
					nmyzsum[mysel]++;
					nzw[hz][w[m][n]]++;
					nzwsum[hz]++;
				} // for n
			} // for m
				// estimate hyperparameters
			estAlpha();
		} // for iter
	} // for run

	// //////////////// query kernel //////////////////
	// Gibbs kernel
	public void runq(int niterq) {
		// iteration loop
		for (int iter = 0; iter < niterq; iter++) {
			System.out.println(iter);
			// major loop, sequence [m][n]
			for (int m = 0; m < Mq; m++) {
				// component selectors
				int mxsel = -1;
				int mxjsel = -1;
				int mysel = -1;
				int myjsel = -1;
				// minor loop, sequence [m][n]
				for (int n = 0; n < wq[m].length; n++) {
					double psum;
					double u;
					// assign topics
					int hx = xq[m][n];
					int hy = yq[m][n];
					int hz = zq[m][n];
					// decrement counts
					nmxq[m][hx]--;
					mxsel = K * m + hx;
					nmxyq[mxsel][hy]--;
					nmxysumq[mxsel]--;
					mysel = L * m + hy;
					nmyzq[mysel][hz]--;
					nmyzsumq[mysel]--;

					// compute weights
					/*
					 * &p(x_{m,n} \eq x, y_{m,n} \eq y, z_{m,n} \eq z\; |\;\vec
					 * x_{-m,n}, \vec y_{-m,n}, \vec z_{-m,n}, \vec w, \cdot)
					 * \notag\\ &\qquad\propto (n^{-mn}_{m,x} + \alpha ) \cdot
					 * \frac{n^{-mn}_{mx,y} + \gamma_{x, y} }{\sum_{y}
					 * n^{-mn}_{mx,y} + \gamma_{x, y}} \cdot
					 * \frac{n^{-mn}_{my,z} + \delta_{y, z} }{\sum_{z}
					 * n^{-mn}_{my,z} + \delta_{y, z}} \cdot \phi_{z,w}
					 */
					psum = 0;
					// hidden edge
					for (hx = 0; hx < K; hx++) {
						// hidden edge
						for (hy = 0; hy < L; hy++) {
							// hidden edge
							for (hz = 0; hz < J; hz++) {
								mxsel = K * m + hx;
								mxjsel = hx;
								mysel = L * m + hy;
								myjsel = hy;

								pp[hx][hy][hz] = (nmxq[m][hx] + alpha)
										* (nmxyq[mxsel][hy] + gamma[mxjsel][hy])
										/ (nmxysumq[mxsel] + gammasum[mxjsel])
										* (nmyzq[mysel][hz] + delta[myjsel][hz])
										/ (nmyzsumq[mysel] + deltasum[myjsel])
										* phi[hz][wq[m][n]];
								psum += pp[hx][hy][hz];
							} // for h
						} // for h
					} // for h

					// sample topics
					u = rand.nextDouble() * psum;
					psum = 0;
					SAMPLED:
					// each edge value
					for (hx = 0; hx < K; hx++) {
						// each edge value
						for (hy = 0; hy < L; hy++) {
							// each edge value
							for (hz = 0; hz < J; hz++) {
								psum += pp[hx][hy][hz];
								if (u <= psum)
									break SAMPLED;
							} // h
						} // h
					} // h

					// assign topics
					xq[m][n] = hx;
					yq[m][n] = hy;
					zq[m][n] = hz;
					// increment counts
					nmxq[m][hx]++;
					mxsel = K * m + hx;
					nmxyq[mxsel][hy]++;
					nmxysumq[mxsel]++;
					mysel = L * m + hy;
					nmyzq[mysel][hz]++;
					nmyzsumq[mysel]++;

				} // for n
			} // for m
		} // for iter
	} // for runq

	// //////////////// hyperparameters //////////////////
	// update hyperparameters
	public void estAlpha() {
		if (iter < 15) {
			return;
		}
		// component selectors
		int mxsel = -1;
		int mxjsel = -1;
		int mysel = -1;
		int myjsel = -1;
		// Note: assuming non-informative gamma priors (1,0)
		// hyperparameter for theta
		int[] nmxsum = new int[M];
		// all components
		for (int m = 0; m < M; m++) {
			nmxsum[m] = w[m].length;
		} // for m
		double xalpha = DirichletEstimation.estimateAlphaMap(nmx, nmxsum,
				alpha, 1., 0.);
		if (alpha < 2.) {
			alpha = (alpha + xalpha) / 2;
		} // < 2
			// hyperparameter for thetax
			// filter nkt and nktsum though jSel index space.
		int[] mx2j = new int[Mq * K];
		// for parent values
		for (int m = 0; m < M; m++) {
			// for parent values
			for (int hx = 0; hx < K; hx++) {
				mxsel = K * m + hx;
				mxjsel = hx;
				mx2j[mxsel] = mxjsel;
			} // for h
		} // for h
		double[][] xgamma = DirichletEstimation.estimateAlphaMapSub(nmxy,
				nmxysum, mx2j, gamma, 1., 0.);
		if (gamma[0][0] < 2.) {

			// all component groups
			for (int j = 0; j < K; j++) {
				// all values
				for (int t = 0; t < L; t++) {
					gamma[j][t] = (xgamma[j][t] + gamma[j][t]) / 2;
				} // for t
			} // for j
		} // < 2
			// hyperparameter for thetay
		// filter nkt and nktsum though jSel index space.
		int[] my2j = new int[Mq * L];
		// for parent values
		for (int m = 0; m < M; m++) {
			// for parent values
			for (int hy = 0; hy < L; hy++) {
				mysel = L * m + hy;
				myjsel = hy;
				my2j[mysel] = myjsel;
			} // for h
		} // for h
		double[][] xdelta = DirichletEstimation.estimateAlphaMapSub(nmyz,
				nmyzsum, my2j, delta, 1., 0.);
		if (delta[0][0] < 2.) {

			// all component groups
			for (int j = 0; j < L; j++) {
				// all values
				for (int t = 0; t < J; t++) {
					delta[j][t] = (xdelta[j][t] + delta[j][t]) / 2;
				} // for t
			} // for j
		} // < 2
			// hyperparameter for phi
		double xbeta = DirichletEstimation.estimateAlphaMap(nzw, nzwsum, beta,
				1., 0.);
		if (beta < 2.) {
			beta = (beta + xbeta) / 2;
		} // < 2
	} // updateHyper

	// //////////////// perplexity //////////////////
	// calculate perplexity value
	public double ppx() {
		double loglik = 0;
		// component selectors
		int mxsel = -1;
		int mxjsel = -1;
		int mysel = -1;
		int myjsel = -1;
		// compute sequence node parameters
		double[][] thetaq = new double[Mq][K];
		// for parent values
		for (int m = 0; m < Mq; m++) {
			int nmxmq = 0;
			for (int x = 0; x < K; x++) {
				nmxmq += nmxq[m][x];
			} // for x
			for (int x = 0; x < K; x++) {
				thetaq[m][x] = (nmxq[m][x] + alpha) / (nmxmq + alphasum);
			} // for x
		} // for parent h
		double[][] thetaxq = new double[Mq * K][L];
		// for parent values
		for (int m = 0; m < Mq; m++) {
			// for parent values
			for (int hx = 0; hx < K; hx++) {
				for (int y = 0; y < L; y++) {
					mxsel = K * m + hx;
					mxjsel = hx;
					thetaxq[mxsel][y] = (nmxyq[mxsel][y] + gamma[mxjsel][y])
							/ (nmxysumq[mxsel] + gammasum[mxjsel]);
				} // for y
			} // for parent h
		} // for parent h
		double[][] thetayq = new double[Mq * L][J];
		// for parent values
		for (int m = 0; m < Mq; m++) {
			// for parent values
			for (int hy = 0; hy < L; hy++) {
				for (int z = 0; z < J; z++) {
					mysel = L * m + hy;
					myjsel = hy;
					thetayq[mysel][z] = (nmyzq[mysel][z] + delta[myjsel][z])
							/ (nmyzsumq[mysel] + deltasum[myjsel]);
				} // for z
			} // for parent h
		} // for parent h
			// compute ppx
		// major loop, sequence [m][n]
		for (int m = 0; m < Mq; m++) {
			// minor loop, sequence [m][n]
			for (int n = 0; n < wq[m].length; n++) {
				double sum = 0;
				// hidden edge
				for (int hx = 0; hx < K; hx++) {
					// hidden edge
					for (int hy = 0; hy < L; hy++) {
						// hidden edge
						for (int hz = 0; hz < J; hz++) {
							mxsel = K * m + hx;
							mysel = L * m + hy;
							sum += thetaq[m][hx] * thetaxq[mxsel][hy]
									* thetayq[mysel][hz] * phi[hz][wq[m][n]];
						} // for h
					} // for h
				} // for h
				loglik += Math.log(sum);
			} // for n
		} // for m
		return Math.exp(-loglik / Wq);
	} // ppx

	// //////////////// monitor string //////////////////
	// describe class and parameters
	public String toString() {
		return "PAM5:\nm >> (theta[m] | alpha) >> x[m][n]\n\t"
				+ "m, x[m][n] >> (thetax[m,x] | gamma[x]) >> y[m][n]\n\t"
				+ "m, y[m][n] >> (thetay[m,y] | delta[y]) >> z[m][n]\n\t"
				+ "z[m][n] >> (phi[z] | beta) >> w[m][n]\n\t"
				+ "\n\t"
				+ "\n"
				+ String.format(
						"Pam5GibbsSampler: \n"
								+ "M = %d  Mq = %d  W = %d  Wq = %d  \n"
								+ "K = %d  L = %d  J = %d  V = %d  \n"
								+ "alpha = %2.5f  gamma[0][0] = %2.5f  delta[0][0] = %2.5f  beta = %2.5f  ",
						M, Mq, W, Wq, K, L, J, V, alpha, gamma[0][0],
						delta[0][0], beta);
	}

} // Pam5GibbsSampler
