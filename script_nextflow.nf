/*
Ce script nextflow se base sur le fichier suivant : "pipeline_STR.sh" disponible sous ce lien : "https://github.com/PYB-SU/STR/blob/main/pipeline_STR.sh"
Ce script est divisée en plusieurs parties : CALL TRIM map zip, OVERLAP, READS, charONT, SEQUENCE ET METHYLATION
Call trim map zip se base sur le fichier suivant "call_trim_zip_map.sh" disponible sur le github
Overlap se base sur le fichier "extract_overlap_GENE.job"
Reads se base sur le fichier "find_best_overlaping_primer.job"
CharONT se base sur le fichier "charONT.job"
Sequence se base sur le fichier "repeat_sequence_from_each_read.job"
Methylation se base sur le fichier "extract_reads_for_methylation_modkit.job" 
Ces fichiers sont tous disponibles sur le github

Ce script nextflow prend en entrée un dossier "PIPELINE" qui contient différents sous-dossiers : 
"CHARONT" qui contiendra tous les fichiers nécessaires pour la partie CharONT
"${params.GENE}" selon le gène choisis qui contiendra tous les fichiers d'entrée nécessaires pour la partie READS 
"FASTA" contient tous les fichiers fasta nécessaires pour la partie Call trim zip map plus particulèrement bonito qui n'est pas présent sur cette version
"FASTQ" contient tous les fichiers fastq nécessaires pour la partie Call trim zip map
"OUTDIR" contient tous les fichiers de sorties de toutes les parties de ce script
	Dans ce dossier, contient "OVERLAPS" qui contiendra un dossier par gène dont celui-ci contiendra tous les fichiers de sorties de la partie OVERLAP
"REFS" contient tous les fichiers de références nécessaires 

Dans le dossier PIPELINE, doit contenir ces deux fichiers suivants : "reads_overlaping_repeat_region.job", "charONT.R"

L'utilisateur a le choix de commencer par un fichier fastq ou bam, cela est définit dans le workflow selon le params.FASTQ

*/
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Fonction pour sélectionner le modèle en fonction de l'entrée de l'utilisateur
def selectionnerModele(model) {
    def models = [
        "FAST10": "dna_r10.4_e8.1_fast@v3.4",
        "HAC10": "dna_r10.4_e8.1_hac@v3.4",
        "SUP10": "dna_r10.4_e8.1_sup@v3.4",
        "FAST81": "dna_r9.4.1_e8.1_fast@v3.4",
        "HAC81": "dna_r9.4.1_e8.1_hac@v3.3",
        "SUP81": "dna_r9.4.1_e8.1_sup@v3.3",
        "FAST8": "dna_r9.4.1_e8_fast@v3.4",
        "HAC8": "dna_r9.4.1_e8_hac@v3.3",
        "SUP8": "dna_r9.4.1_e8_sup@v3.3",
        "TELO": "dna_r9.4.1_telomere"
    ]
    return models[model]
}

// 1) Définition des paramètres nécessaires pour ce script ainsi que les variables nécessaires dans les autres fichiers impliqués


params.GENE = ""
params.MODEL = ""
params.PIPELINEDIR = ""
params.FASTA ="${params.PIPELINEDIR}/FASTA"
params.FASTQ = "${params.PIPELINEDIR}/FASTQ"
params.OUTDIR = "${params.PIPELINEDIR}/OUTDIR"
params.REFS = "${params.PIPELINEDIR}/REFS"
params.REFS_STR = "${params.PIPELINEDIR}/DMPK/haplo_left_CTG_repeat_100.fasta" 
params.CHARONT = "${params.PIPELINEDIR}/CHARONT"
params.help = false
params.cpus = 5
params.THREADS = 32
params.MOD_NAME = ""
params.REPEATLENGTH = ""
params.CUTOFF = "0.5"
params.MAX_READS = "0"
params.BESTPRIMER = "100"
params.QUALITY = 20
params.IDENTITY = 90
params.GENOME = 'chm13v2.0.fa' 
params.coords = "chr19:48597244-48610042"
//params.QUAL_RANGE = [10, 20, 30, 50, 60, 80, 100].toList()
params.PRIMERPOS = ""
params.CUTOFF = ""
                                                                                                                                                                                        

MODEL = selectionnerModele(params.MODEL)

// Afficher les différentes informations input, output, et les différents modèles
def usage() { 
    println ("nextflow.nf --GENE --MODEL --PIPELINEDIR --MOD_NAME --REPEATLENGTH --PRIMERPOS --CUTOFF")
    println(" INPUT : HERE is the directory which contains all files fast5")
    println(" outdir is the output : BAM, FASTQ in a specific folder")
    println(" - use fast5 files in HERE/FAST5")
    println(" - run bonito with each using MODEL") 
    println(" - outputs tmp*fastq.gz files to HERE FASTQ")
    println(" - use porechop with each outputs final fastq.gz to HERE FASTQ")
    println(" - use minimap2 to map file")
    println(" - outputs bam files and index to HERE BAM")
    println(" submits one job to slurm for each file")
    println(" Models are:")
    println("  1) FAST10 dna_r10.4_e8.1_fast@v3.4")
    println("  2) HAC10 dna_r10.4_e8.1_hac@v3.4")
    println("  3) SUP10 dna_r10.4_e8.1_sup@v3.4")
    println("  4) FAST81 dna_r9.4.1_e8.1_fast@v3.4")
    println("  5) HAC81 dna_r9.4.1_e8.1_hac@v3.3")
    println("  6) SUP81 dna_r9.4.1_e8.1_sup@v3.3")
    println("  7) FAST8 dna_r9.4.1_e8_fast@v3.4")
    println("  8) HAC8 dna_r9.4.1_e8_hac@v3.3")
    println("  9) SUP8 dna_r9.4.1_e8_sup@v3.3")
    println(" 10) TELO dna_r9.4.1_telomere")
	println("For the sequence part : ")
	println ("\toutput reads supporting each allele, region matching the repeat, quality over the matching (mean/min/max), matching length and and number of repeats (size/REPEATLENGTH)")
	println("\tall reads are in a single .fastq with coment in name HAP1/HAP2/OUTLIER")
	println("\tuses reads in fasta files from charONT overlap_GENE.gz_reads_allele_[1|2].fasta and overlap_GENE.gz_reads_outliers.fastq")
	println("\taligns with flanking sequences found in ${params.PIPELINEDIR}/GENE/flank-SIDE-SIZE_rev.fa")// comprendre ce que ca veut dire
	println("\tloops over flanking length 20/25/30") 
	println("\t(may require some work to compute actual number of repeats)")
	println("\tThen does the same with ALL reads (initial fastq file) - irrespective of use in alleles")
	println ("\tThen runs straglr on overlap_GENE.bam") // a enlever peut etre
	println("For the methylation part : ")
	println("compute methylation using modkit")
	println("works on BAM files called with dorado/bonito")
	println("needs BAM files in METHYL-ORIG")
}

// Vérification des paramètres
if (params.help) {
    usage()
    exit(1)
} else if (!params.GENE) {
    println("Le gène n'est pas spécifié")
        println("get GENE coords by looking in PIPELINEDIR/GENE/GENE.specs")
        println("output file overlap_GENE.bam in HERE/OVERLAPS dir and index")
        println("loops on all BAM files in HERE/BAM")
    usage()
    exit(1)
}else if (!params.MODEL) {
    println("Le modèle à utiliser pour le traitement des fichiers n'est pas renseigné")
    usage()
    exit(1)
} else if (!params.PIPELINEDIR) {
    println("Le chemin absolu du dossier contenant les fichiers fasta/fastq n'est pas fourni")
    usage()
    exit(1)
}  else if (!params.OUTDIR) {
    println("Le chemin absolu du dossier contenant les fichiers de sortie n'est pas fourni")
    usage()
    exit(1)
} else if (!params.REFS) {
    println("Le chemin absolu du dossier contenant les références nécessaires n'est pas fourni")
    usage()
    exit(1)
}else if (!params.PRIMERPOS) {
    println("PRIMERPOS must be provided")
    usage()
    exit(1)
}else if (!params.CUTOFF) {
    println("CUTOFF must be provided")
    usage()
    exit(1)
} 




// Générer un channel à partir des fichiers fastq.gz dans le répertoire spécifié
Channel_fastq = Channel.fromPath("${params.FASTQ}/*.fastq.gz").map { fastq_gz ->
    def id = fastq_gz.getBaseName().split("\\.")[0] // Extraction de l'ID à partir du nom du fichier
    return [id, fastq_gz]
}.view()

// Générer un channel à partir des fichiers bam dans le répertoire spécifié
Channel_bam_files = Channel.fromPath("${params.OUTDIR}/*.bam").map { bam_file ->
    def id = bam_file.getBaseName().split("\\.")[0] // Extraction de l'ID à partir du nom du fichier
    return [id, bam_file]
}.view()

quality = Channel.fromList(["10", "20", "30", "50", "60", "80", "100"])

/////////////////////////////////////////
//CALL trim map zip
////////////////////////////////////////

// Processus porechop
/*
Ce process prend en entrée un ID d'un fichier et le chemin absolu du fichier
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : "${id}.porechopped.fastq.gz"
Les fichiers de sorties sont dans le dossier OUTDIR
=> porechop : pré-traitement des données pour nettoyer les fichiers fastq
*/
process porechop {
    input:
    tuple val(id), path(fastq_gz) 
    
    output:
    tuple val(id), path("${id}.porechopped.fastq.gz")
    
    script:
    """
    porechop -t ${params.THREADS} -i ${fastq_gz} -o ${id}.porechopped.fastq.gz
    cp ${id}.porechopped.fastq.gz ${params.OUTDIR}
    """
}

// Processus mapping
/*
Ce process prend en entrée un ID d'un fichier et le chemin absolu du fichier : "${id}.porechopped.fastq.gz"
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : "${id}.bam"
Les fichiers de sorties sont dans le dossier OUTDIR
=> minimap2 : alignement de séquences d'ADN longues des fichiers fastq
*/
process mapping {
   // publishDir "${OUTDIR}", mode: 'copy'
    
    input:
    tuple val(id), path(fastq_gz)
    
    output:
    tuple val(id), path("${id}.bam")
    
    script:
    """
    minimap2 -a -x map-ont -t ${params.THREADS} ${params.REFS}/chm13v2.0.mmi ${fastq_gz} | samtools sort -@ ${params.THREADS} -o ${id}.bam
    cp ${id}.bam ${params.OUTDIR}
   
    """
}

//////////////////////////////////////////
// OVERLAP
/////////////////////////////////////////
/*
// Processus get_coords
//process qui ne fonctionne pas, j'attend d'avoir plus de détails sur le pipelindir et ce que cette commande est censé faire exactement
process get_coords {

    input: 
    val(gene)
    
    path(coords)
    output:
    path(coords)

//changer le chemin d'accès, rajouter les coords dans un fichier
    script:
    """

    coords=\$(awk '{print \$2":"\$3"-"\$4}' "${params.PIPELINEDIR}/${gene}/${gene}.specs")
    echo \${coords} >> "${params.PIPELINEDIR}/${gene}/coords.txt"
    """
}
*/

// Processus extract_overlap
/*
Ce process prend en entrée un ID d'un fichier et le chemin absolu du fichier : "${id}.bam"
En sortie, e fichier de sortie avec un nouveau renommage : "${id}_overlap_${params.GENE}.bam"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE
=> samtools view est une commande utilisée pour convertir, filtrer et afficher les fichiers SAM/BAM/CRAM, facilitant l'extraction de sous-ensembles de données d'alignement.
samtools index génère un fichier d'index pour le fichier BAM sortie lors du process mapping, permettant un accès rapide et efficace aux données d'alignement pour des positions spécifiques dans le génome.
*/
process extract_overlap {
    input:
    tuple val(id), path(bam_file)

    output:
    path("${id}_overlap_${params.GENE}.bam")
    
    
    script:
    """
    mkdir -p ${params.OUTDIR}/OVERLAPS

    mkdir -p ${params.OUTDIR}/OVERLAPS/${params.GENE}
    
    samtools index ${bam_file}
    samtools view -@ ${params.THREADS} -b ${bam_file} ${params.coords} > ${id}_overlap_${params.GENE}.bam
    
    cp ${id}_overlap_${params.GENE}.bam ${params.OUTDIR}/OVERLAPS/${params.GENE}
    """
}
//samtools index ${id}_overlap_${params.GENE}.bam ${id}_overlap_${params.GENE}.bai
// Process merge_bams
/*
Ce process prend en entrée le chemin absolu du fichier : "${id}_overlap_${params.GENE}.bam"
En sortie, le fichier de sortie avec un nouveau renommage : "overlap_${params.GENE}_orig.bam"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE
=> samtools merge combine plusieurs fichiers BAM sorties du process extract_overlaping en un seul fichier BAM "overlap_${params.GENE}_orig.bam", en conservant l'ordre des lectures.
samtools index génère un fichier d'index pour le fichier BAM "overlap_${params.GENE}_orig.bam", permettant un accès rapide et efficace aux données d'alignement pour des positions spécifiques dans le génome. Ce fichier d'index se trouve dans le work.
*/
process merge_bams {
    input:
    path(bam_files)

    output: 
    path ("overlap_${params.GENE}_orig.bam")

    script:
    """

    samtools merge -@ ${params.THREADS} -o overlap_${params.GENE}_orig.bam ${bam_files}
	cp overlap_${params.GENE}_orig.bam ${params.OUTDIR}/OVERLAPS/${params.GENE}
    samtools index overlap_${params.GENE}_orig.bam
    cp overlap_${params.GENE}_orig.bam ${params.OUTDIR}/OVERLAPS/${params.GENE}
    """
}

// Process extract_reads
/*
Ce process prend en entrée le chemin absolu du fichier : "overlap_${params.GENE}_orig.bam"
En sortie, les fichiers de sorties avec un nouveau renommage : "overlap_${params.GENE}_rev.fastq.gz", "overlap_${params.GENE}_dir.fastq.gz"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE
=> samtools view est une commande utilisée pour convertir, filtrer et afficher les fichiers BAM, facilitant l'extraction de sous-ensembles de données d'alignement.
=> samtools convertit le fichier BAM "overlap_${params.GENE}_orig.bam" en fichiers FASTQ "overlap_${params.GENE}_rev.fastq.gz","overlap_${params.GENE}_dir.fastq.gz", utilisant le nombre spécifié de threads pour accélérer le processus.
=> seqtk réinverser la séquence et la qualité des lectures dans un fichier FASTQ
*/
process extract_reads {
    input:
    path(bam_file)

    output: 
    tuple path("overlap_${params.GENE}_rev.fastq.gz"), path("overlap_${params.GENE}_dir.fastq.gz")

    script:
    """
    samtools view -@ ${params.THREADS} -b -f 16 -F 32 -F 64 -F 128 -F 256 -F 512 -F 1024 ${bam_file} | samtools fastq -@ ${params.THREADS} | bgzip -@ ${params.THREADS} > overlap_${params.GENE}_rev.fastq.gz
    cp overlap_${params.GENE}_rev.fastq.gz ${params.OUTDIR}/OVERLAPS/${params.GENE} 
    samtools view -@ ${params.THREADS} -b -F 16 -F 32 -F 64 -F 128 -F 256 -F 512 -F 1024 ${bam_file} | samtools fastq -@ ${params.THREADS} | seqtk seq -r | bgzip -@ ${params.THREADS} > overlap_${params.GENE}_dir.fastq.gz
    cp overlap_${params.GENE}_dir.fastq.gz ${params.OUTDIR}/OVERLAPS/${params.GENE}
    """
}

// Process combine_directions
/*
Ce process prend en entrée les chemins absolus du fichier : "overlap_${params.GENE}_rev.fastq.gz", "overlap_${params.GENE}_dir.fastq.gz"
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : "overlap_${params.GENE}_Q0.fastq.gz"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE

*/
process combine_directions {
    input:
    tuple path(rev_reads), path(dir_reads)
	
    output: 
        path ("overlap_${params.GENE}_Q0.fastq.gz")

   script:
   """ 
    	zcat ${rev_reads} ${dir_reads} | bgzip > overlap_${params.GENE}_Q0.fastq.gz
        cp overlap_${params.GENE}_Q0.fastq.gz ${params.OUTDIR}/OVERLAPS/${params.GENE}
   """
}

// Process filter_quality
/*
Ce process prend en entrée les chemins absolus du fichier : "overlap_${params.GENE}_rev.fastq.gz", "overlap_${params.GENE}_dir.fastq.gz"
En sortie,le fichier de sortie avec un nouveau renommage : "overlap_${params.GENE}_Q${qual}.fastq.gz"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE
=> pour filtrer les séquences des fichiers FASTQ en entrée "overlap_${params.GENE}_rev.fastq.gz", "overlap_${params.GENE}_dir.fastq.gz" selon une qualité moyenne {params.quality} et sauvegarde les séquences filtrées dans un fichier FASTQ compressé en sortie "overlap_${params.GENE}_Q${params.quality}.fastq.gz
*/
process filter_quality {
    input:
    tuple path(Q0_fastq_gz), val(qual) 
	
    output: 
    //each params.QUAL_RANGE.collect { path ->
    path ("overlap_${params.GENE}_Q${qual}.fastq.gz")
   // }
	//qual_range = Channel.fromList(['10', '20', '30', '50', '60', '80', '100'])
 	//qual_range = ['10', '20', '30', '50', '60', '80', '100']

   script:
   """ 
   bbduk.sh in=${Q0_fastq_gz} out=overlap_${params.GENE}_Q${qual}.fastq.gz maq=${qual}
   cp overlap_${params.GENE}_Q${qual}.fastq.gz ${params.OUTDIR}/OVERLAPS/${params.GENE}
   """
}

// Process realign_reads
/*
Ce process prend en entrée un ID d'un fichier et le chemin absolu du fichier : "overlap_${params.GENE}_Q${qual}.fastq.gz"
En sortie,les fichiers de sorties avec un nouveau renommage : "overlap_${params.GENE}_samtools_Q${params.quality}.bam", "overlap_${params.GENE}_minimap2_Q${params.quality}.bam"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE
=> Minimap2 est utilisé pour aligner les séquences d'ADN contre un génome de référence (chm13v2.0.mmi).
*/
process realign_reads {
    input:
    path(fastq_files)

    output:
    path("overlap_${params.GENE}_samtools_Q${params.quality}.bam")
	//path("overlap_${params.GENE}_minimap2_Q${params.quality}.bam")
	
    script:
    """
    for fastq_file in \$(ls ${fastq_files}); do

        minimap2 -a -x map-ont "${params.REFS}/chm13v2.0.mmi" \${fastq_file} | samtools view -b -F 32 -F 64 -F 128 -F 256 -F 512 -F 1024 -F 2048 | samtools sort -@ ${params.threads} > overlap_${params.GENE}_minimap2_Q${params.quality}.bam
        cp overlap_${params.GENE}_minimap2_Q${params.quality}.bam ${params.OUTDIR}/OVERLAPS/${params.GENE}
        samtools index overlap_${params.GENE}_minimap2_Q${params.quality}.bam > overlap_${params.GENE}_samtools_Q${params.quality}.bam
        cp overlap_${params.GENE}_samtools_Q${params.quality}.bam ${params.OUTDIR}/OVERLAPS/${params.GENE}
    done
    """
}

/////////////////////////////////////////////
//READS						REFAIRE URGENCE, cette partie prend un fichier en entrée : "reads_overlapping_repeat_region.job"
/////////////////////////////////////////////


// Process find_best_overlaping_reads
/*
Ce process prend en entrée un ID d'un fichier et le chemin absolu du fichier : "overlap_${params.GENE}_Q0.fastq.gz"
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : "find_best_overlaping_reads.log", "overlap_${GENE}_Q0.$direction.sam" ce fichier est le fichier de sortie du fichier "reads_overlapping_repeat_region.job"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE
=> reads_overlapping_repeat_region.job : fichier pris en entrée et est utilisé durant ce process
*/
process find_best_overlaping_reads {
    input:
	path(bam_file)
	val GENE
	val PRIMERPOS
	val PIPELINEDIR
	
    output:
    file "find_best_overlaping_reads.log"

    script:
    """
    LOG="find_best_overlaping_reads.log"
    echo "looking for best PRIMERPOS" > \${LOG}
    echo "PRIMERPOS NB_READS" >> \${LOG}

    
    if [ ! -f "${bam_file}" ]; then
        echo "file ${bam_file} does not exist" >> \${LOG}
        echo "check output of previous" >> \${LOG}
        exit 1
    fi

    MAX_READS=0
    BESTPRIMER=100

    for PRIMERPOS in \$(seq 50 50 1000); do
        NB_READS_PRIMER=\$(sh "${params.PIPELINEDIR}/reads_overlaping_repeat_region.job" "${params.PIPELINEDIR}" "${params.GENE}" "${params.PRIMERPOS}" "${params.CUTOFF}")
        echo "${params.PRIMERPOS} \${NB_READS_PRIMER}" >> \${LOG}
        if [ \${NB_READS_PRIMER} -gt \${MAX_READS} ]; then
            BESTPRIMER=${params.PRIMERPOS}
            MAX_READS=\${NB_READS_PRIMER}
        fi
    done

    echo "best is \${BESTPRIMER}" >> \${LOG}

    zcat "overlap_${params.GENE}_Q0.fastq.gz" > "overlap_${params.GENE}_Q0.fastq"
    cp "overlap_${params.GENE}.fastq" "${params.OUTDIR}/OVERLAPS/${params.GENE}"
    zcat "overlap_${params.GENE}_Q0.fasta.gz" > "overlap_${params.GENE}_Q0.fasta"
    cp "matched_overlap_${params.GENE}_common_reads.ids" "${params.OUTDIR}/OVERLAPS/${params.GENE}"

    echo "${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_Q0.fastq.gz copied to ${params.OUTDIR}/${params.GENE}" >> \${LOG}
    echo "${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_Q0.fasta.gz copied to ${params.OUTDIR}/${params.GENE}" >> \${LOG}

    rm -rf "${params.PIPELINEDIR}/${params.GENE}/tmp_*"

    echo "removed ${params.PIPELINEDIR}/${params.GENE}/tmp_*" >> \${LOG}
    """
}



////////////////////////////////////////
// charONT                                Cette partie prend en entrée un fichier "CharONT.R"
////////////////////////////////////////


// Process charONT
/*
Ce process prend en entrée le chemin absolu du fichier : 
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : "
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/GENE
=> CharONT.R : fichier pris en entrée et est utilisé durant ce process
*/
process charONT {
	input:
	path(charont)
	
	output: 
	
	
	script: 
	"""
    if [ ! -d ${params.OUTDIR}/${params.GENE} ]; then
        mkdir ${params.OUTDIR}/${params.GENE}
	Rscript ${params.CHARONT}/CharONT.R ${params.OUTDIR}/${params.GENE} ${params.GENE}
	"""
}
/*
# echo "copy "${params.CHARONT}/config_CharONT_${params.MODEL}.R" to "${params.OUTDIR}/${params.GENE}/config_CharONT.R" > ${params.LOG}
#cp "${params.CHARONT}/config_CharONT_${params.MODEL}.R ${params.H=}/${params.GENE}/config_CharONT.R"
# echo "removing ${params.HERE}/${params.GENE}/overlap_${params.GENE}.gz" >> ${params.LOG}
# rm -Rf "${params.HERE}/${params.GENE}/overlap_${params.GENE}"
# echo "running charONT on ${params.HERE}/${params.GENE}" >> ${params.LOG}
# cd "${params.HERE}/${params.GENE}"
*/
////////////////////////////////////
// SEQUENCE                           
////////////////////////////////////
/*

process CopyConfig {
    input:
    val(PIPELINEDIR)
    val(GENE)
    path(OUTDIR)

    output:
    file "${params.CHARONT}/config_CharONT.R"

    script:
    """
    cp "${params.PIPELINEDIR}/config_CharONT_${params.MODEL}.R" "${params.OUTDIR}/${params.GENE}/config_CharONT.R"
    """
}
//rm -Rf ${params.OUTDIR}/${params.GENE}/overlap_${params.GENE} VERIFIER CE CHEMIN
process RemoveOverlap {
    input:
    val(OUTDIR)
    val(GENE) 

    script:
    """
    rm -Rf ${params.OUTDIR}/${params.GENE}/overlap_${params.GENE}
    """
}

process RunCharONT {
    input:
    val(PIPELINEDIR)
    val(GENE)
    path(CHARONT)  

// refaire le sc
    script:
    """
    cd ${params.PIPELINEDIR}/${params.GENE}
    Rscript ${params.CHARONT}/CharONT.R ${params.OUTDIR}/${params.GENE} ${params.GENE}
    """
}
*/
// Process SeqtkFasta
/*
Ce process prend en entrée le chemin absolu du fichier : 
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : "overlap_${params.GENE}_reads_outliers.fasta"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE
=> Cette commande utilise Seqtk pour extraire toutes les séquences d'un fichier FASTQ et les convertir en un fichier FASTA. 
*/
process SeqtkFasta {
    input:
    file(fastq) 

    output:
    file "overlap_${params.GENE}_reads_outliers.fasta"

    script:
    """
    seqtk seq -A overlap_${params.GENE}_reads_outliers.fastq > overlap_${params.GENE}_reads_outliers.fasta
    """
}

// Process ProcessFlanking
/*
Ce process prend en entrée le chemin absolu du fichier : "overlap_${params.GENE}_reads_outliers.fasta"
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : "overlap_${params.GENE}_${ALLELE}_${SIZE}_proper.fasta"
Les fichiers de sorties sont dans le dossier GENE qui se trouve dans PIPELINE/OUTDIR/OVERLAPS/GENE
=> Cette commande utilise Seqtk pour extraire toutes les séquences d'un fichier FASTQ et les convertir en un fichier FASTA. 
*/
process ProcessFlanking {
    input:
    val(SIZE)
    val(ALLELE)
    file(file_fasta) 
 
    output:
    file "overlap_${params.GENE}_${ALLELE}_${SIZE}_proper.fasta"

    script:
    """
    for SIDE in left right; do
        water -gapopen 10 -gapextend 1 -outfile stdout ${params.STR_PIPELINEDIR}/${params.GENE}/flank-\${SIDE}-${SIZE}_rev.fa overlap_${params.GENE}_reads_${ALLELE}.fasta -aformat simple -stdout | tr '/' ' ' | awk -v side=\${SIDE} '{if (\$0=="") {} else if (\$2=="Score:") {score=\$3} else if (\$2=="Length:") {lenlocalg=\$3} else if (\$2=="Identity:") {identity=\$3} else if (\$2=="Gaps:") {gaps=\$3} else if (\$2=="2:") {read=\$3} else if (substr(\$1,1,12)==substr(read,1,12)) {if (side=="right") {print "@"read" "\$4" "score" "lenlocalg" "identity" "gaps} else if (side=="left") {print "@"read" "\$2" "score" "lenlocalg" "identity" "gaps}}}' | grep -v -e "^@ " | sort -k1 > overlap_${params.GENE}_${ALLELE}_\${SIDE}${SIZE}.pos
    done
    join overlap_${params.GENE}_${ALLELE}_right${SIZE}.pos overlap_${params.GENE}_${ALLELE}_left${SIZE}.pos > overlap_${params.GENE}_${ALLELE}_${SIZE}.pos
    awk -v size="${SIZE}" -v rl="${params.REPEATLENGTH}" -v allele="${ALLELE}" 'BEGIN{for(n=21;n<256;n++)ord[sprintf("%c",n)]=n-33}NR==FNR{read[\$1]=\$1;left[\$1]=\$2-size;qleft[\$1]=\$3;lenalgleft[\$1]=\$4;idleft[\$1]=\$5;gapleft[\$1]=\$6;right[\$1]=\$7+size;qright[\$1]=\$8;lenalgright[\$1]=\$9;idright[\$1]=\$10;gapright[\$1]=\$11;next}{if (FNR % 4 == 1) {read2=\$1} else if (FNR % 4 == 2) {seq=substr(\$1,left[read2],right[read2]-left[read2])} else if (FNR % 4 == 0) {qual=substr(\$1,left[read2],right[read2]-left[read2]); meanqual=0; leftqual=0; rightqual=0; minqual=256-21;maxqual=0;for(n=0;n<size;n++){pos=substr(qual,n,1);leftqual=leftqual+exp(-ord[pos]*log(10)/10);pos=substr(qual,length(qual)-n+1,1);rightqual=rightqual+exp(-ord[pos]*log(10)/10);}for(n=size;n<length(qual)-size;n++){pos=substr(qual,n,1);meanqual=meanqual+exp(-ord[pos]*log(10)/10);if(minqual>ord[pos])minqual=ord[pos];if(maxqual<ord[pos])maxqual=ord[pos];} meanqual=-10*log(meanqual/(1+length(qual)-2*size))/log(10);rightqual=-10*log(rightqual/size)/log(10);leftqual=-10*log(leftqual/size)/log(10); if (length(seq)>0) {proper="READ-MATCH"} else {proper="READ-NO-MATCH"}; print read2" ( "proper" allele: "allele" avgQ: "meanqual" minQ: "minqual" maxQ: "maxqual" leftQ: "leftqual" rightQ: "rightqual" L: "length(qual)-2*size" N: "(length(qual)-2*size)/rl" S: "qleft[read2]" "qright[read2]" A: "lenalgleft[read2]" "lenalgright[read2]" Id: "idleft[read2]" "idright[read2]" G: "gapleft[read2]" "gapright[read2]" )\n"seq"\n+\n"qual}}' overlap_${params.GENE}_${ALLELE}_${SIZE}.pos overlap_${params.GENE}_reads_${ALLELE}.fastq > overlap_${params.GENE}_${ALLELE}_${SIZE}_proper.fastq
    seqtk seq -a  overlap_${params.GENE}_${ALLELE}_${SIZE}_proper.fastq > overlap_${params.GENE}_${ALLELE}_${SIZE}_proper.fasta
    """

}

/////////////////////////////////////////
// METHYLATION								CHANGER LES INPUTS ET LES OUTPUTS
/////////////////////////////////////////

// Process InitializeMethylation
/*
Ce process prend en entrée : 
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : 
Les fichiers de sorties sont dans le dossier METHYLATION qui se trouve dans OUTDIR/GENE/METHYLATION
=> remplacer scratch_ONT avec le dossier correspondant
*/
process InitializeMethylation {
	input:
	val(OUTDIR)
	val(GENE)
	
	output:
    path "methylation_modkit_${params.PIPELINEDIR.replace('/', '_').replace('scratch_ONT_', '')}_${params.GENE}_${params.MOD_NAME}.log" 
    path "${params.OUTDIR}/${params.GENE}/METHYLATION" 
    
    script: // CRÉATION D'UN LOG ET D'UN DOSSIER METHYLATION + DIR_METHYL
    """
    
    if [ ! -d ${params.PIPELINEDIR}/${params.GENE}/METHYLATION ]; then
        mkdir ${params.PIPELINEDIR}/${params.GENE}/METHYLATION
    fi

    PIPELINEDIR_NO_DIR=\$(echo ${params.PIPELINEDIR} | tr '/' '_' | sed 's/scratch_ONT_//')
    LOG=${params.OUTDIR}/${params.GENE}/METHYLATION/methylation_modkit_\${OUTDIR}_${params.GENE}_${params.MOD_NAME}.log

    cd ${params.OUTDIR}/${params.GENE}
    echo "cd to ${params.OUTDIR}/${params.GENE}" > \${LOG}

    DIR_METHYL=${params.OUTDIR}/${params.GENE}/METHYLATION
    echo " \${DIR_METHYL} exists" >> \${LOG}

    # Output variables
    echo "LOG=\${LOG}"
    echo "DIR_METHYL=\${DIR_METHYL}"
    """ 
}

// Process CheckBAMFiles
/*
Ce process prend en entrée : 
En sortie, le fichier de sortie avec un nouveau renommage : 
Les fichiers de sorties sont dans le dossier METHYLATION qui se trouve dans OUTDIR/GENE/METHYLATION
=> 
*/
process CheckBAMFiles {
    input:
    val LOG 
    val DIR_METHYL 

    script:
    """
    DIR_METHYL=${params.OUTDIR}/${params.GENE}/METHYLATION
    BAMMETHYL=\$(ls \${DIR_METHYL}/*_${params.MOD_NAME}.bam)
    NBBAMMETHYL=\$(echo \${BAMMETHYL} | wc -l)

    if [ \${NBBAMMETHYL} -eq 0 ]; then
        echo "no BAM files for ${params.MOD_NAME} in \${DIR_METHYL}" >>\${LOG}
        exit 1
    elif [ \${NBBAMMETHYL} -ge 2 ]; then
        echo "found \${NBBAMMETHYL} for ${params.MOD_NAME} in \${DIR_METHYL} with ${params.MOD_NAME}.bam - exiting " >>\${LOG}
        exit 1
    else
        echo "found 1 BAM files for ${params.MOD_NAME} in \${DIR_METHYL}" >>\${LOG}
    fi

    # Output variables
    echo "BAMMETHYL=\${BAMMETHYL}"
    """
}
// Process ExportTags
/*
Ce process prend en entrée : 
En sortie, e fichier de sortie avec un nouveau renommage : 
Les fichiers de sorties sont dans le dossier METHYLATION qui se trouve dans OUTDIR/GENE/METHYLATION
=> 
*/
pr
process ExportTags {
    input:
    val LOG 
    val DIR_METHYL
    val BAMMETHYL 

    script:
    """
    # module load samtools

    FILE_NO_EXT=\$(basename \${BAMMETHYL} .bam)
    echo "extract MM and ML from \${BAMMETHYL}" >> \${LOG}
    samtools view \${BAMMETHYL} | awk '{printf \$1; for (i=12; i<=NF; ++i) { if (\$i ~ "^MM:Z:|^ML:B:") {printf "\\t"\$i}}; printf "\\n"}' | sort > overlap_${params.GENE}_${params.MOD_NAME}.tag
	cp overlap_${params.GENE}_${params.MOD_NAME}.tag \${DIR_METHYL}
    # Output file
    echo "TAGFILE=\${DIR_METHYL}/overlap_${params.GENE}_${params.MOD_NAME}.tag"
    """
}

// Process ExtractIDs
/*
Ce process prend en entrée le chemin absolu des fichiers bam : 
En sortie, l'ID et le fichier de sortie avec un nouveau renommage : overlap_${params.GENE}_${params.MOD_NAME}.ids
Les fichiers de sorties sont dans le dossier METHYLATION qui se trouve dans OUTDIR/GENE/METHYLATION
=> 
*/
process ExtractIDs {
    input:
    val LOG 
    
    output:
    path("overlap_${params.GENE}_${params.MOD_NAME}.ids")
    
    script:
    """
    BAMFILE=\$(ls ${params.OUTDIR}/OVERLAPS/${params.GENE}/*.bam)
    NBBAMFILE=\$(echo \${BAMFILE} | wc -l)

    if [ \${NBBAMFILE} -eq 0 ]; then
        echo "no BAM file in ${params.OUTDIR}/OVERLAPS/${params.GENE}" >> \${LOG}
        exit 1
    else
        echo "working with \${BAMFILE}" >> \${LOG}
    fi

    # Obtain IDs from overlap_GENE file
    samtools view overlap_${params.GENE}.bam | awk '{print \$1}' > overlap_${params.GENE}_${params.MOD_NAME}.ids
	cp overlap_${params.GENE}_${params.MOD_NAME}.ids ${params.OUTDIR}/${params.GENE}/METHYLATION
	
    # Output file
    echo "IDFILE=${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_${params.MOD_NAME}.ids"
    """
}

// Process ExtractSAM
/*
Ce process prend en entrée le chemin absolu des fichiers BAM : 
En sortie, les fichiers de sorties avec un nouveau renommage : "overlap_${params.GENE}_${params.MOD_NAME}.sam", "overlap_${params.GENE}_${params.MOD_NAME}.hdr"
Les fichiers de sorties sont dans le dossier METHYLATION qui se trouve dans OUTDIR/GENE/METHYLATION
=> 
*/
process ExtractSAM {
    input:
    val LOG 
    val IDFILE 

    script:
    """
   # module load samtools


    BAMFILE=\$(ls ${params.OUTDIR}/OVERLAPS/${params.GENE}/*.bam)
    for FILE in \${BAMFILE}; do
        samtools view \${FILE} | fgrep -w -f \${IDFILE} | sort >> overlap_${params.GENE}_${params.MOD_NAME}.sam
        cp overlap_${params.GENE}_${params.MOD_NAME}.sam ${params.OUTDIR}/${params.GENE}/METHYLATION
        samtools view -H \${FILE} > overlap_${params.GENE}_${params.MOD_NAME}.hdr
        cp overlap_${params.GENE}_${params.MOD_NAME}.hdr ${params.OUTDIR}/${params.GENE}/METHYLATION
    done


    """
}

// Process ExtractSAM
/*
Ce process prend en entrée le chemin absolu des fichiers fastq : FASTQFILE à changer 
En sortie, les fichiers de sorties avec un nouveau renommage : "overlap_${params.GENE}_${params.MOD_NAME}.sam", "overlap_${params.GENE}_${params.MOD_NAME}.hdr"
Les fichiers de sorties sont dans le dossier METHYLATION qui se trouve dans OUTDIR/GENE/METHYLATION
=> 
*/
process ConvertToFASTQ {
    input:
    val LOG 
    val SAMFILE
    val HEADERFILE
	
	ouput:
	overlap_${params.GENE}_${params.MOD_NAME}_map_STR.sam
	overlap_${params.GENE}_${params.MOD_NAME}_map_STR.hdr
	
    script:
    """
    # module load samtools

    FASTQFILE=${params.OUTDIR}/${params.GENE}/overlap_${params.GENE}_${params.MOD_NAME}.fastq
    SAMFILE=${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_${params.MOD_NAME}.sam
    HEADERFILE=${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_${params.MOD_NAME}.hdr
    if [ "${params.GENE}" == "DMPK" ]; then
        echo "cat \${HEADERFILE} \${SAMFILE} | samtools sort | samtools fastq > \${FASTQFILE}" >> \${LOG}
        cat \${HEADERFILE} \${SAMFILE} > tmp\${SAMFILE}
        samtools fastq tmp\${SAMFILE} > overlap_${params.GENE}_${params.MOD_NAME}.fastq
        cp overlap_${params.GENE}_${params.MOD_NAME}.fastq ${params.OUTDIR}/${params.GENE}/METHYLATION
        
        #module load minimap2

        minimap2 -a -x map-ont ${params.REF_STR} \${FASTQFILE} | sort > overlap_${params.GENE}_${params.MOD_NAME}_map_STR.sam
        cp overlap_${params.GENE}_${params.MOD_NAME}_map_STR.sam ${params.OUTDIR}/${params.GENE}/METHYLATION
        grep -e "^@SQ" -e "^@PG" ${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_${params.MOD_NAME}_map_STR.sam > overlap_${params.GENE}_${params.MOD_NAME}_map_STR.hdr
        cp overlap_${params.GENE}_${params.MOD_NAME}_map_STR.hdr ${params.OUTDIR}/${params.GENE}/METHYLATION
    fi

    # Output file
    echo "FASTQFILE=\${FASTQFILE}"
    """
}

// Process JoinTagsWithSAM
/*
Ce process prend en entrée le chemin absolu des fichiers fastq : FASTQFILE à changer 
En sortie, les fichiers de sorties avec un nouveau renommage : "overlap_${params.GENE}_${params.MOD_NAME}.bam"
Les fichiers de sorties sont dans le dossier METHYLATION qui se trouve dans OUTDIR/GENE/METHYLATION
=> 
*/
process JoinTagsWithSAM {
    input:
    val LOG 
    val SAMFILE 
    val TAGFILE 
    val HEADERFILE 

    script:
    """
    # module load samtools

    METHYLBAMFILE=${params.OUTDIR}/${params.GENE}/METHYLATION/overlap_${params.GENE}_${params.MOD_NAME}.bam

    if [ -f \${TAGFILE} ]; then
        echo "joining \${SAMFILE} \${TAGFILE} in overlap_${params.GENE}_${params.MOD_NAME}.bam" >> \${LOG}
        join \${SAMFILE} \${TAGFILE} -t \$'\\t' > overlap_${params.GENE}_${params.MOD_NAME}
        cat \${HEADERFILE} tmp_overlap_${params.GENE}_${params.MOD_NAME} | samtools sort > \${METHYLBAMFILE}
        rm -f ${params.OUTDIR}/${params.GENE}/METHYLATION/overlap_${params.GENE}_${params.MOD_NAME}
        samtools index overlap_${params.GENE}_${params.MOD_NAME}.bam
        cp overlap_${params.GENE}_${params.MOD_NAME}.bam ${params.OUTDIR}/${params.GENE}/METHYLATION
    fi

    # Output file
    echo "METHYLBAMFILE=\${METHYLBAMFILE}"
    """
}
*

// Process PileupMethylation
/*
Ce process prend en entrée le chemin absolu des fichiers fastq : "overlap_${params.GENE}_${params.MOD_NAME}.bam"
En sortie, les fichiers de sorties avec un nouveau renommage : "overlap_${params.GENE}_${params.MOD_NAME}_CPGfiltered.bedmethyl"
Les fichiers de sorties sont dans le dossier METHYLATION qui se trouve dans OUTDIR/GENE/METHYLATION
=> 
*/
process PileupMethylation {
    input:
    val LOG 
    val METHYLBAMFILE

    script:
    """
    # module load modkit

    BEDMETHYLFILE=${params.OUTDIR}/${params.GENE}/METHYLATION/overlap_${params.GENE}_${params.MOD_NAME}_CPGfiltered.bedmethyl
    modkit pileup \${METHYLBAMFILE} \${BEDMETHYLFILE} --cpg --combine-strands --ref ${params.REFS}/${params.GENOME} --only-tabs

    BEDMETHYLFILE=${params.OUTDIR}/${params.GENE}/METHYLATION/overlap_${params.GENE}_${params.MOD_NAME}_unfiltered.bedmethyl
    modkit pileup \${METHYLBAMFILE} \${BEDMETHYLFILE} --ref ${params.REFS}/${params.GENOME} --only-tabs

    # Output files
    echo "BEDMETHYLFILE_CPG=\${BEDMETHYLFILE_CPG}"
    echo "BEDMETHYLFILE_UNFILTERED=\${BEDMETHYLFILE_UNFILTERED}"
    """
}

////////////////////////////////////////
// WORKFLOW
////////////////////////////////////////


workflow {

	if ( params.FASTQ ) {
// CALL TRIP ZIP MAP
		Channel_fastq | porechop | mapping
		
// OVERLAP
		coords_path = file("${params.PIPELINEDIR}/${params.GENE}/${params.GENE}.specs")

		//get_coords(params.GENE, coords_path )
    	extract_overlap(mapping.out) 
    	merge_bams(extract_overlap.out.toList())
    	extract_reads(merge_bams.out)
    	//qual_range = Channel.fromList(['10', '20', '30', '50', '60', '80', '100'])
    	combine_directions(extract_reads.out)
    	filter_quality(combine_directions.out.combine(quality).view())

    	realign_reads(filter_quality.out.view())
  
// READS
    	find_best_overlaping_reads(combine_directions.out, params.GENE, params.PRIMERPOS, params.PIPELINEDIR)

/*      
// CHARONT   
    	charONT(params.CHARONT)
    
// SEQUENCE
    	CopyConfig(params.PIPELINEDIR, params.GENE, params.OUTDIR)
    	RemoveOverlap(params.PIPELINEDIR, params.GENE)
    	RunCharONT(params.PIPELINEDIR,params.GENE, params.OUTDIR)
    
    	Channel_file_fastq = Channel.fromPath("${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_reads_outliers.fastq")
    	SeqtkFasta(Channel_file_fastq)
    	channel_size = Channel.from([30, 40, 50, 60, 80, 100])
    	channel_allele = Channel.from(['allele_1', 'allele_2', 'outliers', 'all'])
    	channel_fasta = Channel.fromPath("${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_reads_*.fasta")
    	ProcessFlanking(channel_size, channel_allele,channel_fasta)

// METHYLATION
    	InitializeMethylation(params.PIPELINEDIR, params.GENE)

    	BAMFiles_init_LOG = InitializeMethylation.out.find { it.startsWith('LOG') }.split('=')[1]
    	BAMFiles_init_DIR_METHYL = InitializeMethylation.out.find { it.startsWith('DIR_METHYL') }.split('=')[1]    
    
    	CheckBAMFiles(BAMFiles_init_LOG, BAMFiles_init_DIR_METHYL)
    
    	BAMFiles_LOG = CheckBAMFiles.out.find { it.startsWith('LOG') }.split('=')[1]
    	BAMFiles_DIR_METHYL = CheckBAMFiles.out.find { it.startsWith('DIR_METHYL') }.split('=')[1]
    	BAMFiles_BAMMETHYL = CheckBAMFiles.out.find { it.startsWith('BAMMETHYL') }.split('=')[1]
    
    	ExportTags(BAMFiles_LOG, BAMFiles_DIR_METHYL, BAMFiles_BAMMETHYL)
    
    	Export_LOG = ExportTags.out.find { it.startsWith('LOG') }.split('=')[1]
    
		ExtractIDs(Export_LOG)

		ExtractIDS_LOG = ExtractIDs.out.find { it.startsWith('LOG') }.split('=')[1]
		ExtractsIDS_IDFILE = ExtractIDs.out.find { it.startsWith('IDFILE') }.split('=')[1]
		ExtractSAM (ExtractIDS_LOG, ExtractIDS_IDFILE)

		ExtractSAM_LOG = ExtractSAM.out.find { it.startsWith('LOG') }.split('=')[1]
		ExtractSAM_SAMFILE = ExtractSAM.out.find { it.startsWith('SAMFILE') }.split('=')[1]
		ExtractSAM_HEADERFILE = ExtractSAM.out.find { it.startsWith('HEADERFILE') }.split('=')[1]	
		ConvertToFASTQ(	ExtractSAM_LOG,	ExtractSAM_SAMFILE,	ExtractSAM_HEADERFILE)

		ToFASTQ_LOG = ConvertToFASTQ.out.find { it.startsWith('LOG') }.split('=')[1]
		ToFASTQ_SAMFILE = ConvertToFASTQ.out.find { it.startsWith('LOG') }.split('=')[1]
		ToFASTQ_TAGFILE = ExportTags.out.find { it.startsWith('TAGFILE') }.split('=')[1]
		ToFASTQ_HEADERFILE = ConvertToFASTQ.out.find { it.startsWith('HEADERFILE') }.split('=')[1]
		JoinTagsWithSAM(ToFASTQ_LOG,ToFASTQ_SAMFILE,ToFASTQ_TAGFILE,ToFASTQ_HEADERFILE)
	
		Methyl_LOG = JoinTagsWithSAM.out.find { it.startsWith('LOG') }.split('=')[1]
		Methyl_METHYLBAMFILE = JoinTagsWithSAM.out.find { it.startsWith('METHYLBAMFILE') }.split('=')[1]
		PileupMethylation(Methyl_LOG, Methyl_METHYLBAMFILE)
*/
	} else {
		
 
    	//Channel_fastq | porechop | mapping
// OVERLAP
    	get_coords(params.GENE)
    	extract_overlap(get_coords.out.collect(), Channel_bam_files)
    	merge_bams(extract_overlap.out.collect())
    	extract_reads(merge_bams.out)
    	filter_quality(extract_reads.out)
    	realign_reads(filter_quality.out)
    
// READS
    	find_best_overlapping_reads(params.PIPELINEDIR, params.GENE)
/*    
// CHARONT
    	charONT(params.CHARONT)

//SEQUENCE
    	//CopyConfig(params.PIPELINEDIR, params.GENE, params.OUTDIR)
    	//RemoveOverlap(params.PIPELINEDIR, params.GENE)
    	//RunCharONT(params.PIPELINEDIR,params.GENE, params.OUTDIR)
    
    	Channel_file_fastq = Channel.fromPath("${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_reads_outliers.fastq")
    	SeqtkFasta(Channel_file_fastq)
    	channel_size = Channel.from([30, 40, 50, 60, 80, 100])
    	channel_allele = Channel.from(['allele_1', 'allele_2', 'outliers', 'all'])
    	channel_fasta = Channel.fromPath("${params.PIPELINEDIR}/${params.GENE}/overlap_${params.GENE}_reads_*.fasta")
    	ProcessFlanking(channel_size, channel_allele,channel_fasta)

// METHYLATION
    	InitializeMethylation(params.PIPELINEDIR, params.GENE)

    	BAMFiles_init_LOG = InitializeMethylation.out.find { it.startsWith('LOG') }.split('=')[1]
    	BAMFiles_init_DIR_METHYL = InitializeMethylation.out.find { it.startsWith('DIR_METHYL') }.split('=')[1]    
    
    	CheckBAMFiles(BAMFiles_init_LOG, BAMFiles_init_DIR_METHYL)
    
    	BAMFiles_LOG = CheckBAMFiles.out.find { it.startsWith('LOG') }.split('=')[1]
    	BAMFiles_DIR_METHYL = CheckBAMFiles.out.find { it.startsWith('DIR_METHYL') }.split('=')[1]
    	BAMFiles_BAMMETHYL = CheckBAMFiles.out.find { it.startsWith('BAMMETHYL') }.split('=')[1]
    
    
    	ExportTags(BAMFiles_LOG, BAMFiles_DIR_METHYL, BAMFiles_BAMMETHYL)
    
    	Export_LOG = ExportTags.out.find { it.startsWith('LOG') }.split('=')[1]
    
		ExtractIDs(Export_LOG)

		ExtractIDS_LOG = ExtractIDs.out.find { it.startsWith('LOG') }.split('=')[1]
		ExtractsIDS_IDFILE = ExtractIDs.out.find { it.startsWith('IDFILE') }.split('=')[1]
		ExtractSAM (ExtractIDS_LOG, ExtractIDS_IDFILE)
	
		ExtractSAM_LOG = ExtractSAM.out.find { it.startsWith('LOG') }.split('=')[1]
		ExtractSAM_SAMFILE = ExtractSAM.out.find { it.startsWith('SAMFILE') }.split('=')[1]
		ExtractSAM_HEADERFILE = ExtractSAM.out.find { it.startsWith('HEADERFILE') }.split('=')[1]	
		ConvertToFASTQ(	ExtractSAM_LOG,	ExtractSAM_SAMFILE,	ExtractSAM_HEADERFILE)

		ToFASTQ_LOG = ConvertToFASTQ.out.find { it.startsWith('LOG') }.split('=')[1]
		ToFASTQ_SAMFILE = ConvertToFASTQ.out.find { it.startsWith('LOG') }.split('=')[1]
		ToFASTQ_TAGFILE = ExportTags.out.find { it.startsWith('TAGFILE') }.split('=')[1]
		ToFASTQ_HEADERFILE = ConvertToFASTQ.out.find { it.startsWith('HEADERFILE') }.split('=')[1]
		JoinTagsWithSAM(ToFASTQ_LOG,ToFASTQ_SAMFILE,ToFASTQ_TAGFILE,ToFASTQ_HEADERFILE)
	
		Methyl_LOG = JoinTagsWithSAM.out.find { it.startsWith('LOG') }.split('=')[1]
		Methyl_METHYLBAMFILE = JoinTagsWithSAM.out.find { it.startsWith('METHYLBAMFILE') }.split('=')[1]
		PileupMethylation(Methyl_LOG, Methyl_METHYLBAMFILE)
*/
	}


}

