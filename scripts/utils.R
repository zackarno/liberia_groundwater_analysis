####make good colors##########
colv<-c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#333333ff", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E","#BF5B17", "#333333ff", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6")
colv2<-colv[1:15]

colv3<-c( "#386CB0", "#BEAED4", "#FDC086", "#d9f4c7ff","#E7298A","#6dcc29ff", "#fc913eff", "#666666", "#41841cff", "#3bc8ecff", "#7570B3",   "#fcfe42ff","#BF5B17","#333333ff", "#e70b1eff")
colv3<-as.vector(c( c.acidob="#386CB0", c.actino="#BEAED4", c.bacteroid="#FDC086", c.chloro="#d9f4c7ff",c.creno="#E7298A",c.firm="#6dcc29ff", c.gemm="#fc913eff", c.nc10="#666666", c.nitro="#41841cff", c.od1="#3bc8ecff", c.op3="#7570B3",   c.plancto="#fcfe42ff",c.proteob="#BF5B17",c.syn="#333333ff", c.verru="#e70b1eff"))

one_percent_colors<-c("#52095aff", "#bf87c1ff", "#6aae9dff", "#FFFF99", "#386CB0", "#F0027F", "#f3d0b4ff", "#b4d023ff","#067e06ff","#fead33ff","#dead6eff", "#e4202eff", "#f06c59ff","#2b90aeff","#000000","#000000")

trace_elements<-c("V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Sr", "Mo",
                  "Pb")


colv<-c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#333333ff", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E","#BF5B17", "#333333ff", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6")
colv2<-colv[1:15]
cat(paste(shQuote(colv2, type="cmd"), collapse=", "))

colv3<-c( c.acidob="#386CB0", c.actino="#BEAED4", c.bacteroid="#FDC086", c.chloro="#d9f4c7ff",c.creno="#E7298A",c.firm="#6dcc29ff", c.gemm="#fc913eff", c.nc10="#666666", c.nitro="#41841cff", c.od1="#3bc8ecff", c.op3="#7570B3",   c.plancto="#fcfe42ff",c.proteob="#BF5B17",c.syn="#333333ff", c.verru="#e70b1eff")

#define them independently####
c.acidob="#386CB0"; c.actino="#BEAED4"; c.bacteroid="#FDC086"; c.chloro="#d9f4c7ff";c.creno="#E7298A";c.firm="#6dcc29ff"; c.gemm="#fc913eff"; c.nc10="#666666"; c.nitro="#41841cff"; c.od1="#3bc8ecff"; c.op3="#7570B3";   c.plancto="#fcfe42ff";c.proteob="#BF5B17";c.syn="#333333ff"; c.verru="#e70b1eff"

vcols<-as.vector(c(c.acidob, c.actino,  "#FFFF99",c.bacteroid, c.chloro="#d9f4c7ff", c.creno="#E7298A",  c.cyano="orange", c.elusi="#1B9E77",c.eury= "#D95F02", c.fcp="green", c.firm, c.gemm, c.nc10="darkgrey", c.nitro, c.od1, c.op1="#fc546cff", c.op3, c.plancto, c.proteob, c.syn, c.verru,"#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6"))

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}



convert_major_ions_to_mmol<- function(df){
  df %>%
    mutate(
      Ca_mmol= Ca/40.078,
      Mg_mmol= Mg/24.305,
      Na_mmol= Na/22.9897,
      K_mmol= K/39.0983,
      HCO3_mmol= HCO3/((1.008+12.011+(15.999*3))),
      SO4_mmol= SO4/(32.06+(15.999*4)),
      Cl_mmol= Cl/35.45,
      NO3_mmol= NO3/14.007,
      SiO2_mmol=SiO2/(28.085+(15.999*2))
      ) %>%
    select(ends_with("mmol")) %>%
    mutate(H2O_mmol=(1E6/18.02)-rowSums(.))
}
