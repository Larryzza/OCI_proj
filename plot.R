plot_ct_fit <- function(params_df, global_pars, indiv_data, ctalpha=0.01,
                        specialcol="red", basecol="blue", xlim=c(NA,NA), 
                        vlabel="observation", nvlabel="Non-Variant",ntraces=100){
  with(as.list(global_pars),{
    params_df %>% 
      mutate(id_clean=as.integer(id_clean)) %>%
      group_by(id_clean) %>%
      sample_n(ntraces) %>% 
      ungroup() %>%
      ggplot() + 
      # Plot traces:
      geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
      geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
      geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
      geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
      # Plot data:
      geom_point(data=indiv_data, aes(x=t, y=y,col=factor(special)), size=0.5) + 
      scale_color_manual(values=c("1"=specialcol,"0"=basecol), labels=c("1"=vlabel,"0"=nvlabel)) + 
      theme_minimal() + 
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
            legend.title=element_blank(), 
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), 
            axis.text.y=element_text(size=8)) + 
      labs(x="Time since min Ct (days)", y="Ct") + 
      scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), 
                      sec.axis=sec_axis(~convert_Ct_logGEML(.), 
                                        name=expression(log[10]~RNA~copies/ml))) + 
      scale_x_continuous(limits=xlim) + 
      facet_wrap(~id) # Change to id_clean to obscure identities
  })
}

convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
  out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
  return(out) 
}
