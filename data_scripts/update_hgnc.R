updateHGNCdict = function(){

  custom_url = "https://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_status&col=gd_locus_type&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_enz_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=gd_mgd_id&col=gd_other_ids_list&col=gd_pubmed_ids&col=gd_pub_refseq_ids&col=gd_ccds_ids&col=gd_vega_ids&col=md_eg_id&col=md_mim_id&col=md_refseq_id&col=md_prot_id&col=md_ensembl_id&col=md_vega_id&col=md_ucsc_id&col=md_mgd_id&col=md_rgd_id&col=md_rna_central_ids&col=md_lncipedia&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"

  payload = httr::GET(custom_url)

  txt = httr::content(payload,"text")

  hugo_df = read.table(text = txt, stringsAsFactors=F, quote="", comment.char="", sep="\t",
    header = T)

  #build named vector to act as hash object
  #(actual hash object of size required results in segfault)

  syns = plyr::dlply(hugo_df,plyr::.(Approved.Symbol),function(row){
    keys = c(unlist(strsplit(row$Synonyms,",")),
      unlist(strsplit(row$Previous.Symbols,",")),
      unlist(strsplit(row$Accession.Numbers,",")),
      unlist(strsplit(row$Ensembl.Gene.ID,",")),
      unlist(strsplit(row$UniProt.ID.supplied.by.UniProt.,",")),
      unlist(strsplit(row$UCSC.ID.supplied.by.UCSC.,",")),
      unlist(strsplit(row$RefSeq.IDs,",")))
    keys = gsub(" ","",keys)
    keys = unique(keys)
    return(keys)
  })

  hgnc_df = genesetr::dfcols.tochar(stack(syns))


  #ensure none of the keys (alternate symbols) are
  #identical to the values (approved symbols)
  hgnc_dict = c(hgnc_df$ind,unique(hgnc_df$ind))
  names(hgnc_dict) = c(hgnc_df$values,unique(hgnc_df$ind))

  devtools::use_data(hugo_df, overwrite = T)
  devtools::use_data(hgnc_dict, overwrite = T)

}
