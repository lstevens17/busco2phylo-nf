nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
params.outdir = "busco2phylo-${date}"

params.astralDir = "/software/team301/Astral"

busco_dir = Channel.fromPath( "${params.busco_dir_path}", type: 'dir' )

process busco2fasta {
	publishDir "${params.outdir}", mode: 'copy'

	input:
		path(busco_dir)

	output:
		path("b2f_output/*.faa"), emit: input_fastas

	script:
	"""
		busco2fasta.py -b $busco_dir
	"""

}

process align_fastas {
	publishDir "${params.outdir}/alignments", mode: 'copy'
	
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

	input:
		path input_alignments

	output:
		path("*.treefile"), emit: gene_trees

	script:
		"""
		iqtree -s $input_alignments -nt ${task.cpus}
		"""
}

process run_astral {
	publishDir "${params.outdir}/astral", mode: 'copy'

	input:
		path gene_trees
	
	output:
		path("astral_tree.nwk"), emit: astral_tree

	script:
		"""
		cat ${gene_trees} >all_gene_trees.txt
		java -jar $params.astralDir/astral.5.7.4.jar -i all_gene_trees.txt -o astral_tree.nwk
		"""
	
}

process trim_alignments {
	publishDir "${params.outdir}/trimmed_alignments", mode: 'copy'
	
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
	
	input:
		path trimmed_alignments
	
	output:
		path("supermatrix.fa"), emit: supermatrix
		
	script:
		"""
		catfasta2phyml.pl -c -f ${trimmed_alignments} >supermatrix.fa
		"""
}

process infer_supermatrix_tree {
	publishDir "${params.outdir}/iqtree", mode: 'copy'
	
	input: 
		path supermatrix
	
	output: 
		path("*")
	
	script: 
		"""
		iqtree -s ${supermatrix} -bb 1000 -nt ${task.cpus}
		"""
}

process estimate_astral_BLs {
	publishDir "${params.outdir}/astral", mode: 'copy'
	
	input: 
		tuple path(astral_tree), path(supermatrix)

	output:
		path("*")
	
	script: 
		"""
		iqtree -s ${supermatrix} -nt ${task.cpus} -te ${astral_tree}
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
