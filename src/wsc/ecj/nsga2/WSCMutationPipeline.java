package wsc.ecj.nsga2;

import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Parameter;
import wsc.data.pool.Service;
import wsc.problem.WSCInitializer;

public class WSCMutationPipeline extends BreedingPipeline {

	private static final long serialVersionUID = 1L;

	@Override
	public Parameter defaultBase() {
		return new Parameter("wscmutationpipeline");
	}

	@Override
	public int numSources() {
		return 1;
	}

	@Override
	public int produce(int min, int max, int start, int subpopulation, Individual[] inds, EvolutionState state,
			int thread) {

		int n = sources[0].produce(min, max, start, subpopulation, inds, state, thread);

		// Individual original = (Individual) inds[start];

		if (!(sources[0] instanceof BreedingPipeline)) {
			inds[start] = (Individual) (inds[start].clone());
		}

		if (!(inds[start] instanceof SequenceVectorIndividual))
			// uh oh, wrong kind of individual
			state.output
					.fatal("WSCMutationPipeline didn't get a SequenceVectorIndividual. The offending individual is: "
							+ inds[start]);

		WSCInitializer init = (WSCInitializer) state.initializer;

		// Perform mutation
		SequenceVectorIndividual tree = (SequenceVectorIndividual) inds[start].clone();
		int indexA = init.random.nextInt(tree.genome.length);
		int indexB = init.random.nextInt(tree.genome.length);
		swapServices(tree.genome, indexA, indexB);
		inds[start] = tree;
		inds[start].evaluated = false;

		return n;
	}

	private void swapServices(Service[] genome, int indexA, int indexB) {
		Service temp = genome[indexA];
		genome[indexA] = genome[indexB];
		genome[indexB] = temp;
	}

}
