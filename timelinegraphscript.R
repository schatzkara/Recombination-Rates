years = 1928:2012

penicillin = vector("numeric",85)
tetracycline = vector("numeric",85)
methicillin = vector("numeric",85)
erythromycin = vector("numeric",85)
gentamicin = vector("numeric",85)
ceftazidime = vector("numeric",85)
vancomycin = vector("numeric",85)
levofloxacin = vector("numeric",85)
imipenem = vector("numeric",85)
linezolid = vector("numeric",85)
daptomycin = vector("numeric",85)
ceftaroline = vector("numeric",85)

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1928 && j <= 1940) {
    penicillin[[i]] = 1
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1950 && j <= 1959) {
    tetracycline[[i]] = 2
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1953 && j <= 1968) {
    erythromycin[[i]] = 3
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1960 && j <= 1962) {
    methicillin[[i]] = 4
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1967 && j <= 1979) {
    gentamicin[[i]] = 5
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1972 && j <= 1988) {
    vancomycin[[i]] = 6
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1985 && j <= 1998) {
    imipenem[[i]] = 7
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1985 && j <= 1987) {
    ceftazidime[[i]] = 8
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 1996 && j <= 1997) {
    levofloxacin[[i]] = 9
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 2000 && j <= 2001) {
    linezolid[[i]] = 10
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 2003 && j <= 2005) {
    daptomycin[[i]] = 11
  }
  j = j+1
  i = i+1
}

j = 1928
i = 1
while (i <= 85) {
  if (j >= 2010 && j <= 2011) {
    ceftaroline[[i]] = 12
  }
  j = j+1
  i = i+1
}

dataMat = matrix(c(penicillin,tetracycline,erythromycin,methicillin,gentamicin,vancomycin,imipenem,ceftazidime,levofloxacin,linezolid,daptomycin,ceftaroline),nrow=85)
colnames(dataMat) = c("penicillin","tetracycline","erythromycin","methicillin","gentamicin","vancomycin","imipenem","cefatzidime","levofloxacin","linezolid","daptomycin","ceftaroline")
rownames(dataMat) = years
dataMat = data.frame(dataMat)

i <- 1
j <- 1
while (j <= 12) {
  while (i <= 85) {
    if (dataMat[i,j] == 0) {
      dataMat[i,j] <- NA
    } 
    i = i+1  
  }
  j = j+1
  i = 1
}

plot = ggplot(data = dataMat,aes(x=years))
#plot = plot+geom_line(aes(x = years, y = dataMat$penicillin,color="penicillin"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$tetracycline,color="tetracycline"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$erythromycin,color="erythromycin"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$methicillin,color="methicillin"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$gentamicin,color="gentamicin"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$vancomycin,color="vancomycin"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$imipenem,color="imipenem"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$cefatzidime,color="cefatzidime"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$levofloxacin,color="levofloxacin"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$linezolid,color="linezolid"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$daptomycin,color="daptomycin"),size=3,show.legend=F)
plot = plot+geom_line(aes(x = years, y = dataMat$ceftaroline,color="ceftaroline"),size=3,show.legend=F)
plot = plot+theme_classic()
plot = plot+theme(axis.line.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
plot = plot+labs(x = NULL,color = "Antibiotic")
plot = plot+scale_x_continuous(breaks = seq(from = 1950, to = 2010, by = 10),limits = c(1950,2015))
plot = plot+theme(axis.text.x = element_text(size=24))

print(plot)
  