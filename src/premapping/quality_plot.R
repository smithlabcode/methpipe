quality_plot <- function(datafile, baseplot = "", qualplot = "",
						make_baseplot = TRUE, make_qualplot = FALSE)
{
	cols = c("red", "green", "blue", "yellow", "black");
	cols = rainbow(5);

	d = read.table(datafile, header=T);
	read_len = length(d[[1]]);	
	total  = sum(d[1,6:10]);
	
	for (i in 6:10) {d[i] = d[[i]] / total;}

	for (i in 7:10) {d[i] = d[[i]] + d[[i-1]];}

	if (make_baseplot || baseplot != "") 
	{
		plot(c(0,0), col="white", 
				xlim=c(0, read_len), ylim=c(0,1),
				xlab="bases", ylab="frequency",
				main="Base composition along reads");

		for (i in 1:read_len)
		{
			rect(i-1.5, 0, i - 0.5, d[i, 6], col=cols[1]);
			for (j in 6:9) 
				rect(i-1.5, d[i, j], i - 0.5, d[i, j+1], col=cols[j-4]);
		}

		if (baseplot != "")
			dev.copy2eps(file=baseplot);
	}

	if (make_qualplot || qualplot != "")
	{
		plot(d[,5], main="Quality score along reads", 
			xlab="bases", ylab="quality score")

		if (qualplot != "")
			dev.copy2eps(file=qualplot);	
	}	
}

if (length(commandArgs) == 3)
	quality_plot(datafile = commandArgs[1],
				 baseplot = commandArgs[2],
				 qualplot = commandArgs[3]);
	
