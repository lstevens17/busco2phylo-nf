nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "busco2phylo-${date}"

busco_dir = Channel.fromPath( "${params.busco_dir_path}", type: 'dir' )

process busco2fasta {
        publishDir "${params.outdir}", mode: 'copy'
        queue 'normal'

        input:
                path(busco_dir)

        output:
                path("b2f_output/*.faa"), emit: input_fastas

        script:
        """
                python3 /lustre/scratch116/tol/teams/team301/projects/lepidoptera_genomics/tree_building/busco2fasta/busco2fasta.py -b $busco_dir
        """

}

process align_fastas {
        publishDir "${params.outdir}/alignments", mode: 'copy'
        conda = "/nfs/users/nfs_l/ls30/miniconda3/envs/orthology_env/"

        queue 'small'

        input:
                path(input_fastas)

        output:
                path("*.aln"), emit: input_alignments

        script:
        """
                cut -f1 -d'.' ${input_fastas} >${input_fastas}.reformatted
                mafft --auto ${input_fastas}.reformatted > ${input_fastas}.reformatted.aln
        """
}

process infer_gene_trees {
	publishDir "${params.outdir}/gene_trees", mode: 'copy'
	conda = "/nfs/users/nfs_l/ls30/miniconda3/envs/orthology_env/"
	queue 'normal'
	cpus 2

	input:
                path input_alignments

        output:
                path("*.treefile"), emit: gene_trees

        script:
        """
              iqtree -s $input_alignments -nt 2
        """
}

process run_astral {
	publishDir "${params.outdir}/astral", mode: 'copy'
	queue 'normal'

	input:
		path gene_trees
	
	output:
		path("astral_tree.nwk"), emit: astral_tree

	script:
	"""
		cat ${gene_trees} >all_gene_trees.txt
		java -jar /lustre/scratch116/tol/teams/team301/projects/lepidoptera_genomics/tree_building/Astral/astral.5.7.4.jar -i all_gene_trees.txt -o astral_tree.nwk
	"""
	
}

process trim_alignments {
	publishDir "${params.outdir}/trimmed_alignments", mode: 'copy'
	conda = "/nfs/users/nfs_l/ls30/miniconda3/envs/orthology_env/"
	queue 'small'
	
	input: 
		path input_alignments
	
	output:
		path("*.trimal"), optional: true, emit: trimmed_alignments
	
	script:
	"""
		trimal -in ${input_alignments} -out ${input_alignments}.trimal -gt 0.8 -st 0.001 -resoverlap 0.75 -seqoverlap 80
	"""
}

process concatenate_alignments {
	publishDir "${params.outdir}/supermatrix", mode: 'copy'
	queue 'normal'
	
	input:
                path trimmed_alignments

        output:
                path("supermatrix.fa"), emit: supermatrix

        script:
        """
		/lustre/scratch116/tol/teams/team301/projects/lepidoptera_genomics/tree_building/catfasta2phyml.pl -c -f ${trimmed_alignments} >supermatrix.fa
        """
}

process infer_supermatrix_tree {
	publishDir "${params.outdir}/iqtree", mode: 'copy'
	conda = "/nfs/users/nfs_l/ls30/miniconda3/envs/orthology_env/"
	queue 'long'
	cpus 16 
	
	input: 
		path supermatrix
	
	output: 
		path("*")
	
	script: 
	"""
		iqtree -s ${supermatrix} -bb 1000 -nt 16
	"""
}

process estimate_astral_BLs {
	publishDir "${params.outdir}/astral", mode: 'copy'
	conda = "/nfs/users/nfs_l/ls30/miniconda3/envs/orthology_env/"
	queue 'long'
	cpus 16
	
	input: 
		tuple path(astral_tree), path(supermatrix)

	output:
		path("*")
	
	script: 
	"""
		iqtree -s ${supermatrix} -nt 16 -te ${astral_tree}
	"""
}

workflow {

        busco2fasta(busco_dir) 
        busco2fasta.out.flatten() | align_fastas | infer_gene_trees
	align_fastas.out | trim_alignments
	trim_alignments.out.collect() | concatenate_alignments | infer_supermatrix_tree
	infer_gene_trees.out.collect() | run_astral
	run_astral.out.combine(concatenate_alignments.out) | estimate_astral_BLs
}
