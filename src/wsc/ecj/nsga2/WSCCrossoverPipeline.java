package wsc.ecj.nsga2;

import java.util.HashSet;
import java.util.Set;

import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Parameter;
import wsc.data.pool.Service;
import wsc.problem.WSCInitializer;

public class WSCCrossoverPipeline extends BreedingPipeline {

	private static final long serialVersionUID = 1L;

	@Override
	public Parameter defaultBase() {
		return new Parameter("wsccrossoverpipeline");
	}

	@Override
	public int numSources() {
		return 2;
	}

	@Override
	public int produce(int min, int max, int start, int subpopulation,
			Individual[] inds, EvolutionState state, int thread) {

		WSCInitializer init = (WSCInitializer) state.initializer;

		Individual[] inds1 = new Individual[inds.length];
		Individual[] inds2 = new Individual[inds.length];

		int n1 = sources[0].produce(min, max, 0, subpopulation, inds1, state, thread);
		int n2 = sources[1].produce(min, max, 0, subpopulation, inds2, state, thread);

        if (!(sources[0] instanceof BreedingPipeline)) {
            for(int q=0;q<n1;q++)
                inds1[q] = (Individual)(inds1[q].clone());
        }

        if (!(sources[1] instanceof BreedingPipeline)) {
            for(int q=0;q<n2;q++)
                inds2[q] = (Individual)(inds2[q].clone());
        }

        if (!(inds1[0] instanceof SequenceVectorIndividual))
            // uh oh, wrong kind of individual
            state.output.fatal("WSCCrossoverPipeline didn't get a SequenceVectorIndividual. The offending individual is in subpopulation "
            + subpopulation + " and it's:" + inds1[0]);

        if (!(inds2[0] instanceof SequenceVectorIndividual))
            // uh oh, wrong kind of individual
            state.output.fatal("WSCCrossoverPipeline didn't get a SequenceVectorIndividual. The offending individual is in subpopulation "
            + subpopulation + " and it's:" + inds2[0]);

        int nMin = Math.min(n1, n2);

        // Perform crossover
        for(int q=start,x=0; q < nMin + start; q++,x++) {
    		SequenceVectorIndividual t1 = ((SequenceVectorIndividual)inds1[x]);
    		SequenceVectorIndividual t2 = ((SequenceVectorIndividual)inds2[x]);

    		// Select two random index numbers as the boundaries for the crossover section
    		int indexA = WSCInitializer.random.nextInt(t1.genome.length);
    		int indexB = WSCInitializer.random.nextInt(t1.genome.length);

    		// Make sure they are different
    		while (indexA == indexB)
    			indexB = WSCInitializer.random.nextInt(t1.genome.length);

    		// Determine which boundary they are
    		int minBoundary = Math.min(indexA, indexB);
    		int maxBoundary = Math.max(indexA, indexB);

    		// Create new genomes
    		Service[] newGenome1 = new Service[t1.genome.length];
    		Service[] newGenome2 = new Service[t2.genome.length];

    		// Swap crossover sections between candidates, keeping track of which services are in each section
    		Set<Service> newSection1 = new HashSet<Service>();
    		Set<Service> newSection2 = new HashSet<Service>();

    		for (int index = minBoundary; index <= maxBoundary; index++) {
    			// Copy section from parent 1 to genome 2
    			newGenome2[index] = t1.genome[index];
    			newSection2.add(t1.genome[index]);

    			// Copy section from parent 2 to genome 1
    			newGenome1[index] = t2.genome[index];
    			newSection1.add(t2.genome[index]);
    		}

    		// Now fill the remainder of the new genomes, making sure not to duplicate any services
    		fillNewGenome(t2, newGenome2, newSection2, minBoundary, maxBoundary);
    		fillNewGenome(t1, newGenome1, newSection1, minBoundary, maxBoundary);

    		// Replace the old genomes with the new ones
    		t1.genome = newGenome1;
    		t2.genome = newGenome2;

	        inds[q] = t1;
	        inds[q].evaluated=false;

	        if (q+1 < inds.length) {
	        	inds[q+1] = t2;
	        	inds[q+1].evaluated=false;
	        }
        }
        return n1;
	}

	private void fillNewGenome(SequenceVectorIndividual parent, Service[] newGenome, Set<Service> newSection, int minBoundary, int maxBoundary) {
		int genomeIndex = getInitialIndex(minBoundary, maxBoundary);

		for (int i = 0; i < parent.genome.length; i++) {
			if (genomeIndex >= newGenome.length + 1)
				break;
			// Check that the service has not been already included, and add it
			if (!newSection.contains(parent.genome[i])) {
				newGenome[genomeIndex] = parent.genome[i];
				// Increment genome index
				genomeIndex = incrementIndex(genomeIndex, minBoundary, maxBoundary);
			}
		}
	}

	private int getInitialIndex(int minBoundary, int maxBoundary) {
		if (minBoundary == 0)
			return maxBoundary + 1;
		else
			return 0;
	}

	private int incrementIndex(int currentIndex, int minBoundary, int maxBoundary) {
		if (currentIndex + 1 >= minBoundary && currentIndex + 1 <= maxBoundary)
			return maxBoundary + 1;
		else
			return currentIndex + 1;
	}
}
