all:
	# to pdf
	R --slave -e 'rmarkdown::render("bioutils.Rmd", output_format ="rmarkdown::pdf_document")'
	# general html documnet --> better scrollable tables
	#R --slave -e 'rmarkdown::render("bioutils.Rmd", output_format ="rmarkdown::html_document")'
	# html document specifically formatted for vignette
	R --slave -e 'rmarkdown::render("bioutils.Rmd", output_format ="rmarkdown::html_vignette")'
	# to github formatted markdown
	R --slave -e 'rmarkdown::render("bioutils.Rmd", output_format ="rmarkdown::github_document")'

clean:
	rm -f bioutils.md bioutils.html bioutils.pdf
	rm -rf bioutils_files/
