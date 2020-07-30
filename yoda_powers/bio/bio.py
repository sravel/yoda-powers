def concat_fasta_files(path_directory):
    """
    Return a fasta dictionnary of concatenation fasta file's find in directory ("fasta", "fa", "fas")

    Warning:
        Sequence on fasta must have the same name

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        path_directory (str): a path to fasta file directory

    Returns:
        :class:`dict`: python dict with the concatenation of fasta filename in path_directory ( file with extention "fa", "fasta", "fas" )

    Raises:
         ValueError: If `path_directory` does not exist.
         ValueError: If `path_directory` is not a valid directory.

    Examples:
        >>> dico_concat = concat_fasta_files('path/to/directory/')
        >>> print(dico_concat)
            {"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
    """
    from pathlib import Path
    from Bio import SeqIO

    if not Path(path_directory).exists():
        raise ValueError(f'ERROR: directory "{path_directory}" does not exist')
    elif not Path(path_directory).is_dir():
        raise ValueError(f'ERROR: "{path_directory} " is not a valid directory')

    path_directory = Path(path_directory).resolve()
    fa_files = path_directory.glob("*.fa")
    fasta_files = path_directory.glob("*.fasta")
    fas_files = path_directory.glob("*.fas")

    output_dico_seqs = {}

    for fasta_file in (*fa_files, *fasta_files, *fas_files):
        # print(fasta_file)
        with open(fasta_file.as_posix(), "rU") as handle:
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        for seq_name in sorted(record_dict.keys()):
            if seq_name not in output_dico_seqs.keys():
                name = record_dict[seq_name].id
                seq = record_dict[seq_name].seq
                output_dico_seqs[name] = seq
            else:
                name = record_dict[seq_name].id
                seq = record_dict[seq_name].seq
                output_dico_seqs[name] += seq

    return output_dico_seqs


def convert_fasta_2_nexus(path_directory, path_directory_out):
    """
    Return the number of fasta file's convert  find in directory ("fasta", "fa", "fas") where are converted

    Warning:
        Sequence on fasta must align and have the same length

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        path_directory (str): a path to fasta file directory
        path_directory_out (str): a directory path to write nexus file

    Returns:
        :class:`int`: the number of file converted

    Raises:
         ValueError: If `path_directory` does not exist.
         ValueError: If `path_directory` is not a valid directory.
         ValueError: If fasta is not align.

    Examples:
        >>> nb_file = convert_fasta_2_nexus('path/to/directory/',' path/to/directory/')
        >>> print(nb_file)
            "4172"
    """
    from pathlib import Path
    from Bio import AlignIO
    from Bio.Nexus import Nexus

    if not Path(path_directory).exists():
        raise ValueError(f'ERROR: directory "{path_directory}" does not exist')
    elif not Path(path_directory).is_dir():
        raise ValueError(f'ERROR: "{path_directory} " is not a valid directory')

    path_directory = Path(path_directory).resolve()
    path_directory_out = Path(path_directory_out).resolve()
    if not path_directory_out.exists():
        path_directory_out.mkdir()

    fa_files = path_directory.glob("*.fa")
    fasta_files = path_directory.glob("*.fasta")
    fas_files = path_directory.glob("*.fas")

    # general variables:
    minimal_record = f'#NEXUS\nbegin data; dimensions ntax=0 nchar=0; format datatype=dna; end;'
    count_convert = 0

    for fasta_file in (*fa_files, *fasta_files, *fas_files):
        try:
            count_convert += 1
            print(Path(fasta_file))
            basename = Path(fasta_file).stem
            alignment = AlignIO.read(fasta_file.as_posix(), format='fasta')
            n = Nexus.Nexus(minimal_record)
            n.alphabet = alignment._alphabet
            for record in alignment:
                n.add_sequence(record.id.replace("-", "_"), str(record.seq))
            n.write_nexus_data(f"{path_directory_out.as_posix()}/{basename}.nex", interleave=False)
        except ValueError as e:
            raise ValueError(f"ERROR on file {fasta_file}, with message: {e}, please check")

    return count_convert


def dict_2_fasta(dico, fasta_out):
    """
    Function that takes a dictionary where key are ID and value Seq, and write a fasta file.

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        dico (dict): python dict with ID in key and Seq on values
        fasta_out (str): a directory path to write nexus file

    Returns:
        :class:`str`: the output fasta file name

    Examples:
        >>> dico = {"Seq1":"ATGCTGCAGTAG","Seq2":"ATGCCGATCGATG","Seq3":"ATGCTCAGTCAGTAG"}
        >>> dict_2_fasta(dico)
        >Seq1
        ATGCTGCAGTAG
        >Seq2
        ATGCCGATCGATG
        >Seq3
        ATGCTCAGTCAGTAG
    """
    from pathlib import Path
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    with open(Path(fasta_out).resolve().as_posix(), "w") as output_handle:
        for seq_name in sorted(dico.keys()):
            seq_obj = Seq(dico[seq_name])
            record = SeqRecord(seq_obj, id=seq_name, name=seq_name, description="")
            SeqIO.write(record, output_handle, "fasta")
    return output_handle.name


def extract_seq_from_fasta(fasta_file, wanted_file, include=True):
    """
    Function to extract sequence from fasta file

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        fasta_file (str): a path to fasta file directory
        wanted_file (str): a path file with id (one per line)
        include (bool, optional): if True keep id on wanted file, else keep id not in wanted file

    Returns:
        :class:`dict`: the fasta dict with sequence extract

    Raises:
         ValueError: If `wanted_file` or `fasta_file` does not exist.
         ValueError: If `wanted_file` or `fasta_file`` is not a valid file.
         ValueError: If `include` is not valid boolean.

    Example:
        >>> dict_sequences = extract_seq_from_fasta(fasta_file, wanted_file)
        >>> dict_sequences
        {'Seq2': SeqRecord(seq=Seq('ATGCCGATCGATG', SingleLetterAlphabet()), id='Seq2', name='Seq2', description='Seq2'
        , dbxrefs=[]), 'Seq3': SeqRecord(seq=Seq('ATGCTCAGTCAGTAG', SingleLetterAlphabet()), id='Seq3', name='Seq3',
        description='Seq3', dbxrefs=[])}
    """
    from Bio import SeqIO
    from pathlib import Path

    wanted_file = Path(wanted_file).resolve()
    fasta_file = Path(fasta_file).resolve()

    if not isinstance(bool(include), bool):
        raise ValueError(f'ERROR: "{include}" must be a boolean value')
    include = bool(include)

    if not wanted_file.exists():
        raise ValueError(f'ERROR: file "{wanted_file}" does not exist')
    elif not wanted_file.is_file():
        raise ValueError(f'ERROR: "{wanted_file} " is not a valid file')

    if not fasta_file.exists():
        raise ValueError(f'ERROR: file "{fasta_file}" does not exist')
    elif not fasta_file.is_file():
        raise ValueError(f'ERROR: "{fasta_file} " is not a valid file')

    with open(wanted_file) as f:
        wanted = set([line.strip() for line in f if line != ""])

    with open(fasta_file, "r") as handle:
        fasta_sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        fasta_sequences_del = fasta_sequences.copy()
    for seq in fasta_sequences:
        if include and seq not in wanted:
            del fasta_sequences_del[seq]
        elif not include and seq in wanted:
            del fasta_sequences_del[seq]
    return fasta_sequences_del


def fasta_2_dict(fasta_file):
    """
    Function that take a file name (fasta), and return a dictionnary of sequence

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        fasta_file (str): a path to fasta file directory

    Returns:
        :class:`dict`: the fasta dict with sequence extract

    Raises:
         ValueError: If `fasta_file` does not exist.
         ValueError: If `fasta_file` is not a valid file.

    Example:
        >>> filename = "sequence.fasta"
        >>> fasta_2_dict(filename)
        {'Seq1': SeqRecord(seq=Seq('ATGCTGCAGTAG', SingleLetterAlphabet()), id='Seq1', name='Seq1', description='Seq1', dbxrefs=[]),
        'Seq2': SeqRecord(seq=Seq('ATGCCGATCGATG', SingleLetterAlphabet()), id='Seq2', name='Seq2', description='Seq2', dbxrefs=[]),
        'Seq3': SeqRecord(seq=Seq('ATGCTCAGTCAGTAG', SingleLetterAlphabet()), id='Seq3', name='Seq3', description='Seq3', dbxrefs=[])}
    """
    from Bio import SeqIO
    from pathlib import Path

    fasta_file = Path(fasta_file).resolve()

    if not fasta_file.exists():
        raise ValueError(f'ERROR: file "{fasta_file}" does not exist')
    elif not fasta_file.is_file():
        raise ValueError(f'ERROR: "{fasta_file} " is not a valid file')

    with open(fasta_file, "r") as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))


def len_seq_2_dict(fasta_file):
    """
    Function that take a file name (fasta), and return a dictionnary with length of sequence

    Notes:
        function need modules:

        - pathlib
        - BioPython

    Arguments:
        fasta_file (str): a path to fasta file directory

    Returns:
        :class:`dict`: the fasta dict with sequence length

    Raises:
         ValueError: If `fasta_file` does not exist.
         ValueError: If `fasta_file`` is not a valid file.

    Example:
        >>> filename = "sequence.fasta"
        >>> len_seq_2_dict(filename)
        {'Seq1': 12, 'Seq2': 13, 'Seq3': 15}
    """
    from Bio import SeqIO
    from pathlib import Path

    fasta_file = Path(fasta_file).resolve()
    dico_lenght = {}

    if not fasta_file.exists():
        raise ValueError(f'ERROR: file "{fasta_file}" does not exist')
    elif not fasta_file.is_file():
        raise ValueError(f'ERROR: "{fasta_file} " is not a valid file')

    with open(fasta_file, "r") as handle:
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    for gene in sorted(record_dict.keys()):
        if record_dict[gene].id not in dico_lenght:
            lenseq = len(record_dict[gene].seq)
            dico_lenght[gene] = int(lenseq)
    return dico_lenght


def nb_seq_files_2_dict(path_directory):
    """
    Function  that take a Path Directory and returna dictionnary with number of sequences in fasta file's

    :param pathDirectory: a directory Path
    :type pathDirectory: Path
    :rtype: dict1(), dict2()
    :return: - contient le nombre de sequences dans les fichiers (key = nom de fichier value =  nombre de sequences)\n
             - contient le nombre de fichier qui ont x sequences (key = nombre de sequence =  nombre de fichier)
    :raise print: print("ERROR: Sequence: "+nameFichier+" allready read") with nameFichier is the current file read.

    Example:
        >>> dico1,dico2 = nbSeqInFile2dict(path/to/directory/)
        >>> print(dict2txt(dico1))
        ./out/gemo10_4497_ortho_rename_add.fasta	58
        ./out/gemo10_6825_ortho_rename_add.fasta	59
        ./out/gemo10_3497_ortho_rename_add.fasta	59
        ./out/gemo10_6254_ortho_rename_add.fasta	59
        >>> print(dict2txt(dico2))
        58	1
        59	3
    """

    dicoNbSeqInFiles = {}
    dicoNbFilesNbSouche = {}
    # récupération des fichiers du repertoire
    if pathDirectory[-1] != "*":
        pathDirectory += "*"
    directoryFiles = glob.glob(pathDirectory)
    # Ouverture des fichiers du repertoire pour stocker les sequences en mémoire
    for fichier in directoryFiles:
        try:
            nameFichier = fichier.split("/")[-1].split(".")[0]
            extentionFichier = fichier.split("/")[-1].split(".")[1]
        except:
            extentionFichier = "directory"
        if extentionFichier in ["fasta", "fa", "fas"]:
            # Ouverture des sequences fasta et chargement dans dictionnaire
            handle = open(fichier, "r")
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            handle.close()
            nbseq = len(record_dict.keys())
            # print(nameFichier+"\t"+str(nbseq))
            if nameFichier not in dicoNbSeqInFiles.keys():
                dicoNbSeqInFiles[nameFichier] = nbseq
            else:
                print("ERROR: Sequence: " + nameFichier + " allready read")

            dicoNbFilesNbSouche[nbseq] = (dicoNbFilesNbSouche.get(nbseq, 1)) + 1
    return dicoNbSeqInFiles, dicoNbFilesNbSouche


class ParseGFF:
    """
    Parser of GFF3 file write in python.
    return an object iterable contain GFFRecord()

    line in GFF3 return:

    Example:
        >>> scaffold_44     prediction      gene    46      6942    0       +       .       ID=gene_1;Name=jgi.p|Mycfi2|180833;portal_id=Mycfi2;proteinId=180833;transcriptId=180833
        >>> GFFRecord(seqid='scaffold_44', source='prediction', type='gene', start=46, end=6942, score=0.0, strand='+', phase=None,
        >>> 		attributes={'portal_id': 'Mycfi2', 'transcriptId': '180833', 'proteinId': '180833', 'Name': 'jgi.p|Mycfi2|180833', 'ID': 'gene_1'}, seq=None, len=6896)

    GFFRecord has attributes can acces with record.value (ex: record.seqid):

    ==========  ===========================================================
    attribute    infos
    ==========  ===========================================================
    seqid       first column of gff3
    source      second column of gff3
    type        third column of gff3 contain type
    start       start position
    end         end position
    score       score
    strand      DNA brin
    phase       phase
    attributes  dict() with key corresponding to GFFAttributes
    seq         if fasta load can add sequence but by default = None
    len         size of sequence
    ==========  ===========================================================

    Example:
        >>> objGFF = ParseGFF(gffFile)
        >>> for record in objGFF.parseGFF3():
        >>> 	print(record.seqid)
        >>> 	if record.type == "mRNA" :
        >>> 		transcriptID = record.attributes["transcriptId"]

    """

    def __init__(self, filename):
        from collections import namedtuple
        # Initialized GeneInfo named tuple. Note: namedtuple is immutable
        self.filename = filename
        self.gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes",
                              "seq", "len"]
        self.GFFRecord = namedtuple("GFFRecord", self.gffInfoFields)

    @staticmethod
    def parseGFFAttributes(self, attributeString):
        """Parse the GFF3 attribute column and return a dict"""
        import urllib
        if attributeString == ".": return {}
        ret = {}
        for attribute in attributeString.split(";"):
            key, value = attribute.split("=")
            ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
        return ret

    def parseGFF3(self):
        """
        A minimalistic GFF3 format parser.
        Yields objects that contain info about a single GFF3 feature.

        Supports transparent gzip decompression.
        """
        import gzip
        import urllib
        # Parse with transparent decompression
        open_fn = gzip.open if self.filename.endswith(".gz") else open
        with open_fn(self.filename) as in_file:
            for line in in_file:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                # If this fails, the file format is not standard-compatible
                assert len(parts) == len(self.gffInfoFields) - 2
                # Normalize data
                normalizedInfo = {
                        "seqid"     : None if parts[0] == "." else urllib.parse.unquote(parts[0]),
                        "source"    : None if parts[1] == "." else urllib.parse.unquote(parts[1]),
                        "type"      : None if parts[2] == "." else urllib.parse.unquote(parts[2]),
                        "start"     : None if parts[3] == "." else int(parts[3]),
                        "end"       : None if parts[4] == "." else int(parts[4]),
                        "len"       : None if parts[4] == "." and parts[3] == "." else int(parts[4]) - int(parts[3]),
                        "score"     : None if parts[5] == "." else float(parts[5]),
                        "strand"    : None if parts[6] == "." else urllib.parse.unquote(parts[6]),
                        "phase"     : None if parts[7] == "." else urllib.parse.unquote(parts[7]),
                        "seq"       : None,
                        "attributes": self.parseGFFAttributes(parts[8])
                }
                # Alternatively, you can emit the dictionary here, if you need mutability:
                #	yield normalizedInfo
                yield self.GFFRecord(**normalizedInfo)
