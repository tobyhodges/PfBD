# Trimming by window method with minimum length cutoff

def trimByWindowWithMinLength(sequence, qualityThreshold, minimumLength):
	for pos in xrange(len(sequence), 1, -5):
		if pos < minimumLength:
			return None
			break
		qScoreIndex = pos-1
		window = sequence.letter_annotations['phred_quality'][qScoreIndex-5:qScoreIndex]
		meanWindowQual = sum(window)/5.0
		if meanWindowQual >= qualityThreshold:
			for windowPos in xrange(pos,pos-5,-1):
				windowQScoreIndex = windowPos-1
				if sequence.letter_annotations['phred_quality'][windowQScoreIndex] >= qualityThreshold:
					trimmedSeq = sequence[:windowPos]
					return trimmedSeq
					break
			break


# Reverse Complement Anti-Sense Sequences

def revCompAntiSense(sequences):
	processed = []
	for sequence in sequences:
		if not str(sequence.seq).startswith('ATG'):
			rcSeq = sequence.reverse_complement(id=sequence.id + '_revcomp')
			processed.append(rcSeq)
		else:
			processed.append(sequence)
	return processed

# Translate Sequences

def translateRecords(sequences):
	translated = []
	for sequence in sequences:
		translatedRec = SeqIO.SeqRecord(seq=sequence.seq.translate(), id=sequence.id + '_translated')
		translated.append(translatedRec)
	return translated