#! /path/to/Rscript --vanilla --default-packages=utils


#Notes: GRanges object imported from different format should be specified as format parameter
# when exported to ensure correct format i.e. format=...

# group is a List of factor, List accessed by: [[i]]; factor should be as.character() to a vector first

listOverlaps <- function(query, subject, minoverlap=0.90) {
	hits <- list();
	qRange <- ranges(query);
	sRange <- ranges(subject);
	for (i in 1:length(qRange@width)) {
		for (j in 1:length(sRange@width)) {
			qstart <- qRange@start[i];
			qoverlap <- as.integer(qRange@width[i] * minoverlap);
			qend <- qstart + qRange@width[i] - 1;
			sstart <- sRange@start[j];
			soverlap <- as.integer(sRange@width[i] * minoverlap);
                        send <- sstart + sRange@width[j] - 1;
			if (sstart <= qstart) {
				if (send >= qend || (send - qstart) >= qoverlap || (send - qstart) >= soverlap) {
					hits <- c(hits, list(c(i,j)));
				}
			} else if (send <= qend || (qend - sstart) >= qoverlap || (qend - sstart) >= soverlap) {
                                        hits <- c(hits, list(c(i,j)));
                        }
		}
	}
	return (hits);
}


annotTransfer <- function(target_file, source_file, minoverlap = 0.90, features = c("gene"), replace = TRUE, runTrinotate = TRUE, seqFile = "../genome_S288C_R64.fsa", outFile = "sample") {

library("GenomicRanges", lib.loc="../lib/GenomicRanges");

library("rtracklayer", lib.loc="../lib/GenomicRanges");

library("seqinr", lib.loc="../lib");

minoverlap <- 0.90

target <- import.gff(target_file, format="gff");

source <- import.gff(source_file, format="gff");


targetData <- as.character(target@elementMetadata@listData$group);
  
targetSeqNames <- as.vector(seqnames(target));
  
targetFeatures <- as.vector(target@elementMetadata@listData$type);


sourceData <- as.character(source@elementMetadata@listData$group);

sourceSeqNames <- as.vector(seqnames(source));

seq <- read.fasta(seqFile, forceDNAtolower = FALSE);
#str(seq);
tempSeq <- list();
cat(sprintf("Running annotation transfer\n"));

sourceFeatures <- as.vector(source@elementMetadata@listData$type);

listLines <- list();

listFeatures <- list();
#name - position mapping used for appending trinotate result
nameMap <- list();
if (replace) {
  trans.dupLog <- list();
  dup.dupLog <- list();
  show(system.time(overlap <- listOverlaps(target, source, minoverlap)));
  q();
  dupCheck <- 0;
  lastAnnotated <- 0;
  annotated <- FALSE;
  output <- GRanges();
  for (o in overlap) {
	 i <- o[[1]];
	 j <- o[[2]];
	 # if a line does not overlap anything, or no annotation has been transferred, add to trinotate list
	 if (i != lastAnnotated) {
		for (k in lastAnnotated + 1:i - 1) {
		  cat(sprintf("No existing annotation found at line %d, added to trinotate list\n", k));
		  #Get ID using regex substring
		  ID <- sub(".*ID=([^;]*).*", "\\1", targetData[k]);
		  cat(sprintf("ID:%s\n", ID));
		  output <- c(output, target[k]);
		  targetTemp <- target[k];
		  targetTemp@elementMetadata@listData$type <- factor(c("exon"));
		  output <- c(output, targetTemp);
		  a <- as.SeqFastadna(get(targetSeqNames[k],seq)[start:end], name = ID, Annot = paste(">", ID));
		  #str(a);
		  if (length(tempSeq) == 0) {
			 tempSeq <- list(a);
		  } else {
			 tempSeq <- c(tempSeq, list(a));
		  }
		  names(tempSeq) <- c(names(tempSeq)[1:length(names(tempSeq)) - 1], ID);
		  nameMap[[ID]] <- length(output@ranges@width);
		}
	 } 
	 if (targetSeqNames[i] == sourceSeqNames[j] && sourceFeatures[j] %in% features) {
		lastAnnotated <- i;
		if (j %in% listLines) {
		  cat(sprintf("Line %d dropped from source: %s\n", j, sourceSeqNames[j]));
		  next;
		}
		listLines <- c(listLines, j);
		source.id <- strsplit(strsplit(sourceData[j], ";")[[1]][1], "=")[[1]][2];
		if (i == dupCheck) {
		  dup.dupLog <- c(dup.dupLog, list(source.id));
		}
		# target transfer not duplicated, check strand
		else {
		  dupCheck <- i;
		  trans.dupLog <- c(trans.dupLog, list(source.id));
		  if (as.vector(target@strand)[i] == as.vector(source@strand)[j]) {
			 output <- c(output, source[j]);
			 name <- sub(".*Name=([^;]*).*", "\\1",sourceData[j]);
			 if (name != sourceData[j]) {
				k = j + 1;
				cat(sprintf("Name: %s\nParent: %s\n", name, sub(".*Parent=([^;]*).*", "\\1",sourceData[k])));
				while ((parent <- sub(".*Parent=([^;]*).*", "\\1",sourceData[k])) != sourceData[k] && parent == name) {
				  cat(sprintf("Transferred line %d from source: child of %s\n", k, parent));
				  output <- c(output, source[k]);
				  listLines <- c(listLines, k);
				  k <- k + 1;
				}
			 }
			 
			 if ((feature <- sub(".*Name=([^;]*).*", "\\1",sourceData[j])) != sourceData[j]) {
				
				listFeatures <- c(listFeatures, list(feature));
			 }
			 cat(sprintf("Transferred annotation from line %d to line %d\n", j, i));
		  } else {
			 cat(sprintf("Line %d marked as ncRNA from line %d\n", i, j));
			 targetTemp <- target[i];
			 targetTemp@elementMetadata@listData$type <- factor(c("ncRNA"));
			 output <- c(output, targetTemp);
			 targetTemp@elementMetadata@listData$type <- factor(c("exon"));
			 output <- c(output, targetTemp);
		  }
		}
	 }
  }
  target <- output;
  targetData <- as.character(target@elementMetadata@listData$group);
  targetSeqNames <- as.vector(seqnames(target));
  targetFeatures <- as.vector(target@elementMetadata@listData$type);
  if (length(dup.dupLog) != 0) {
	 n <- max(length(dup.dupLog), length(trans.dupLog));
	 length(dup.dupLog) <- n;
	 length(trans.dupLog) <- n;
	 cat(sprintf("Writing duplicate log\n"));
	 #str(cbind(as.vector(dup.dupLog),as.vector(trans.dupLog)));
	 #str(cbind(unlist(dup.dupLog), unlist(trans.dupLog)));
	 write.table(cbind(dup.dupLog, trans.dupLog), file = "duplicates.log", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
  
  
} else {
  for (i in 1:length(target@ranges@width)) {
	 cat(sprintf("Processing line: %d\n", i));
	 start <- start(ranges(target))[i];
	 width <- width(ranges(target))[i];
	 end <- end(ranges(target))[i];
	 query <- IRanges(c(start), c(end));
	 overlap <- as.integer(width * minoverlap);
	 #cat(sprintf("Required overlap: %d\n", overlap));
	 #cat(sprintf("Start:%d, End:%d\n", start, end));
	 hits <- findOverlaps(query, ranges(source), minoverlap=overlap);
	 annotated <- FALSE;
	 #show(hits);
	 for (j in subjectHits(hits)) {
		# transfer annotation
		if (targetSeqNames[i] == sourceSeqNames[j] && sourceFeatures[j] %in% features) {
		  if (targetSeqNames[i] %in% listLines) {
			 
			 next;
		  }
		  annotated <- TRUE;
		  listLines <- c(listLines, targetSeqNames[i]);
		  if (as.vector(target@strand)[i] == as.vector(source@strand)[j]) {
			 targetData[i] <- paste( sourceData[j], "(", start(ranges(source))[j], "--", end(ranges(source))[j], ")",";", targetData[i]);
			 cat(sprintf("Transferred annotation from line %d to line %d\n", j, i));
		  } else {
			 # cat(sprintf("Line %d marked as ncRNA\n", i));
			 targetFeatures[i] = "ncRNA";
		  }
		}
	 }
	 if (!annotated) {
		cat(sprintf("No existing annotation found at line %d, added to trinotate list\n", i));
		#Get ID using regex substring
		ID <- sub(".*ID=(.*);.*", "\\1", targetData[i]);
		#str(ID);
		
		a <- as.SeqFastadna(get(targetSeqNames[i],seq)[start:end], name = ID, Annot = paste(">", ID));
		#str(a);
		if (length(tempSeq) == 0) {
		  tempSeq <- list(a);
		} else {
		  tempSeq <- c(tempSeq, list(a));
		}
		names(tempSeq) <- c(names(tempSeq)[1:length(names(tempSeq)) - 1], ID);
		nameMap[[ID]] <- i;
		#cat(sprintf("length: %d\n", length()))
		#str(tempSeq);
	 }
	 
  }
}
#str(tempSeq);


#write.fasta(tempSeq, names = names(tempSeq), file.out = "temp.fasta");

# integrate trinotate result to unannotated regions

#system("./trinotAdd.sh temp.fasta Trinotate.sqlite");
trinotate <- read.delim("trinotate_annotation_report.xls", row.names = NULL);

transnames <- as.vector(trinotate[[2]]);
for (i in 1:length(transnames)) {
	for (j in 4:length(trinotate) - 1) {
		a <- transnames[i];
	#	cat(sprintf("transnames: %s\n", a));
		if (as.vector(trinotate[[j]])[i] != ".") {
			targetData[nameMap[[transnames[i]]]] <- paste(targetData[nameMap[[transnames[i]]]], ";", names(trinotate)[j + 1], ":", as.vector(trinotate[[j]])[i]);
	#	Writing trinotate data for the duplicated "exon" line
			targetData[nameMap[[transnames[i]]] - 1] <- paste(targetData[nameMap[[transnames[i]]]], ";", names(trinotate)[j + 1], ":", as.vector(trinotate[[j]])[i]);
		}
	}
}
target@elementMetadata@listData$group <- factor(targetData);
target@elementMetadata@listData$type <- factor(targetFeatures);
cat(sprintf("Writing output..\n"));
write.table(as.data.frame(target), paste(outFile, ".txt", sep = ""), sep="\t");
write.csv(as.data.frame(target), paste(outFile, ".csv", sep = ""));
export.gff(target, paste(outFile, ".gff", sep = ""));

}

annotTransfer("../tri6.pasa_assemblies.gff3", "../tri6.gene_structures_post_PASA_updates.65674.gff3");
