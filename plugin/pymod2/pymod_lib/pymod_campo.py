from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import numpy
import pymod_sequence_manipulation as pmsm

class CAMPO:
    """
    A class to analyze a protein multiple alignment using the CAMPO algorithm first described in:
        - Paiardini A1, Bossa F, Pascarella S., Nucleic Acids Res. 2005 Jul 1;33(Web Server issue):W50-5.
          CAMPO, SCR_FIND and CHC_FIND: a suite of web tools for computational structural biology.
    """

    def __init__(self,fasta_file_full_path, mutational_matrix = "blosum62", gap_score=-1, gap_gap_score=0, toss_gaps=True):
        """
        Takes as an argument the file path of the .fasta format alignment file.
        """
        self.fasta_file_full_path = fasta_file_full_path
        self.records = list(SeqIO.parse(self.fasta_file_full_path, "fasta"))
        self.num_seq= len(self.records)

        self.sequence_list=[]
        for seq_record in self.records:
            sequence = str(seq_record.seq)
            self.sequence_list.append(sequence)

        # Test that all sequences have the same length.
        same_length = None
        if len(set([len(s) for s in self.sequence_list])) == 1:
            same_length = True
            self.alignment_length = len(self.sequence_list[0])
        else:
            same_length = False

        if not same_length:
            raise Exception("Le sequenze non sono allineate e non e' possibile calcolare i CAMPO scores.")

        # Prepares the substitution matrix.
        self.mutational_matrix = None
        if mutational_matrix == "blosum62":
            self.mutational_matrix = MatrixInfo.blosum62
        elif mutational_matrix == "blosum90":
            self.mutational_matrix = MatrixInfo.blosum90
        elif mutational_matrix == "blosum80":
            self.mutational_matrix = MatrixInfo.blosum80
        elif mutational_matrix == "blosum50":
            self.mutational_matrix = MatrixInfo.blosum50
        elif mutational_matrix == "blosum45":
            self.mutational_matrix = MatrixInfo.blosum45
        elif mutational_matrix == "pam250":
            self.mutational_matrix = MatrixInfo.pam250
        elif mutational_matrix == "pam120":
            self.mutational_matrix = MatrixInfo.pam120
        elif mutational_matrix == "pam30":
            self.mutational_matrix = MatrixInfo.pam30

        # Completes the "other half" of the Biopython matrix.
        for pair in self.mutational_matrix.keys():
            value = self.mutational_matrix[pair]
            reversed_pair = (pair[1],pair[0])
            self.mutational_matrix.update({reversed_pair: value})

        # Adds values for X.
        for r in ("C","S","T","P","A","G","N","D","E","Q","H","R","K","M","I","L","V","F","Y","W","X"):
            x_value = -1
            self.mutational_matrix.update({(r,"X"): x_value, ("X",r): x_value})

        # Adds items for residue-gaps pairs.
        for r in ("C","S","T","P","A","G","N","D","E","Q","H","R","K","M","I","L","V","F","Y","W","X"):
            self.mutational_matrix.update({(r,"-"): gap_score, ("-",r): gap_score})

        # Adds values for gap-gap pair.
        self.mutational_matrix.update({("-","-"): gap_gap_score})

        # If this variable is set to True, positions with too many gaps will not be assigned with a
        # CAMPO score.
        self.toss_gaps = toss_gaps


    def compute_id_matrix(self):
        """
        Compute the identity matrix, necessary to calculate CAMPO scores.
        This should be called just after an object of this class is built.
        """
        self.id_matrix = []
        for i in xrange (self.num_seq-1):
            self.id_matrix.append([])

        for i in range (0, self.num_seq-1):
            for j in range (i+1, self.num_seq):
                identity = pmsm.compute_sequence_identity(self.sequence_list[i], self.sequence_list[j])
                self.id_matrix[i].append(identity)


    def run_CAMPO(self):
        """
        Actually runs the CAMPO algorithm and stores the conservation scores for each column of the
        multiple alignment.
        """

        #############################
        # Italian code starts here. #
        #############################

        self.matrice_somme=[]
        for i in xrange(0, self.num_seq-1):
            self.matrice_somme.append([])

        # Genera una variabile denominatore che conterra' la sommatoria di tutti i valori (1-matrice[i][indice])
        denominatore= 0.0

        # Questo ciclo annidato risolve il numeratore della frazione contenuta nell'algoritmo di CAMPO
        # Inoltre aggiunge a denominatore il valore di (1-matrice[i][indice])
        for i in range (0, self.num_seq-1):

            # Questo indice servira' per richiamare la giusta % identita' dalla matrice delle identita'
            indice=0
            for j in range (i+1, self.num_seq):

                # Lista[] alla fine del ciclo conterra' in maniera ordinata i valori ottenuti dal confronto
                # del primo amminoacido della sequenza i con il primo della sequenza j, del secondo AA di i
                # col secondo AA di j e cosi' via...

                self.lista=[]

                # Questo ciclo confronta il primo amminoacido della sequenza i con il primo amminoacido della sequenza j
                # Ne calcola in Bscorek(ij) e lo divide per (|(Bscorek(ii)|+|Bscorek(jj)|)*(1/2)
                # Bl[AAi,AAj] restituisce il Bscore dello scambio AAi-->AAj
                # Questo risultato viene moltiplicato per 1-%identita'(ij) e successivamente aggiunto a lista[]
                for AAi, AAj in zip (self.sequence_list[i],self.sequence_list[j]):
                    if AAi != '-' or AAj != '-':
                        numeratore=0
                        blosum_term = (abs(self.get_match_score((AAi,AAi)))+abs(self.get_match_score((AAj,AAj)))) * float(0.5)
                        try:
                            numeratore= self.get_match_score((AAi,AAj)) /  (blosum_term) * (1-self.id_matrix[i][indice])
                        except Exception,e:
                            numeratore= self.get_match_score((AAi,AAj)) /  (1) * (1-self.id_matrix[i][indice])

                        self.lista.append(round(float(numeratore),2))
                    else:
                        numeratore=self.get_match_score((AAi,AAi))
                        self.lista.append(round(float(numeratore),2))

                denominatore= denominatore+(1-self.id_matrix[i][indice])
                self.matrice_somme[i].append(self.lista)
                indice=indice+1

                # Questo controllo serve per evitare che effettui una divisione 0/0 nel caso in cui tutte le sequenze siano identiche.
                # Ok sara' comunque pari a zero per ogni colonna dell'allineamento poiche' tutti i valori al numeratore della sommatoria
                # vengono moltiplicati per 1-matrice[i][indice] che nel caso di sequenze identiche e' pari a zero.

                if denominatore == 0:
                    denominatore= 0.01

        # Alla fine del ciclo il primo elemento contenuto in matrice_somme[] sara la matrice contenente i valori ottenuti
        # dal confronto della prima sequenza con la seconda, poi con la terza e cosi via. Il secondo elemento sara' la
        # matrice che contiene i valori dei confronti tra la seconda sequenza con la terza, poi con la quarta e cosi via...

        # Questo ciclo risolve la sommatoria al numeratore dell'algoritmo di CAMPO

        # Il ciclo piu' esterno (i) scorre le varie colonne da sommare. Quello intermedio (z) seleziona quale delle matrici contenute
        # in matrice_somme sto prendendo in considerazione mentre il ciclo (j) identifica le righe di quella matrice.

        self.lista_somme_colonne=[]
        for i in range (0,len(self.sequence_list[1])):
            somma=0
            for z in range(0, len(self.matrice_somme)):
                for j in range (0, len(self.matrice_somme[z])):
                    somma=somma+self.matrice_somme[z][j][i]
            self.lista_somme_colonne.append(somma)

        # Questo ciclo calcola i punteggi Ok assegnati ad ogni colonna K dell'allineamento
        # valori_Ok e' la lista che conterra', in maniera ordinata, tutti i punteggi assegnati alle varie colonne dell'allineamento

        self.campo_scores=[]

        for somma_colonna in self.lista_somme_colonne:
            valore= (1 / (self.num_seq*(self.num_seq-1)*float(0.5))) * float(somma_colonna) / float(denominatore)
            self.campo_scores.append(valore)

        # The fraction of gaps in a column of a multiple alignment in order for it to not
        # be assigned a normalized CAMPO score.
        self.toss_gap_threshold = 0.25
        self.tossed_alignment_positions = []
        if self.toss_gaps:
            for alignment_position in range(self.alignment_length):
                gap_count = 0
                for seq in self.sequence_list:
                    if seq[alignment_position] == "-":
                        gap_count += 1
                gap_fraction = float(gap_count)/float(self.num_seq)
                if gap_fraction >= self.toss_gap_threshold:
                    self.tossed_alignment_positions.append(alignment_position)
                    self.campo_scores[alignment_position] = None

        # Normalize on the maximum value. Tosses alignment positions with too much gaps.
        massimo = max(self.campo_scores)
        self.normalized_campo_scores = []
        for score in self.campo_scores:
            if score != None:
                self.normalized_campo_scores.append(score/massimo)
            else:
                self.normalized_campo_scores.append(None)

        ###########################
        # Italian code ends here. #
        ###########################


    def get_match_score(self,residues_pair):
        score = None
        try:
            score = self.mutational_matrix[residues_pair]
        except KeyError:
            score = 0
        return score


    def get_campo_items_list(self):
        # Tosses alignment positions with too much gaps.
        filtered_scores = filter(lambda x : x != None, self.normalized_campo_scores)
        clist = numpy.array([round(i, 3) for i in filtered_scores])
        min_campo = min(clist)
        max_campo = max(clist)
        bins = numpy.array(numpy.linspace(min_campo, max_campo, num=10))
        inds = numpy.digitize(clist,bins)
        list_of_campo_items = []
        for campo_score, bin_id in zip(clist, inds):
            list_of_campo_items.append({"campo-score": campo_score, "interval": bin_id})
        # Adds back 'None' values for position tossed out because of their high gap content.
        for tossed_position in self.tossed_alignment_positions:
            list_of_campo_items.insert(tossed_position, {"campo-score": None, "interval": None})
        return list_of_campo_items
