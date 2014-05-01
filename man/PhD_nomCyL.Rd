\name{PhD_nomCyL}
\alias{PhD_nomCyL}
\docType{data}
\title{
Data set from Survey on Human Resources in Science and Technology carried out by
Spanish Statistical Office.
}
\description{
The sample data corresponds 681 answers from PhD holders, 
corresponding to people that in 2006 had a doctoral degree and with their
residence in Castilla-Le\'on region in Spain. The data come from Survey on Careers
of doctorate holders(CDH) carried out by Spanish Statistical Office in 2008.
}
\usage{data(PhD_nomCyL)}
\format{
  This data frame contains 681 observation for the following 7 columns:
  \describe{
    \item{MS}{ Marital Status:(1:M(Married),2:MLR(Living in a marriage-like relationship),
				3:SD (Separated or Divorced),4:SW(Widowed or Single)}  
    \item{SECT}{Sector of employment(1:BES(Business Enterprise Sector),
				2:GS (Government Sector),3:HES(Higher Education Sector),
				4:PNP(Private Non Profit))}
    \item{MIN}{Minimum education level required for the principal job:
				(1:mPD(Postdoc),2:mARQ(Advanced Research Qualification),
				3:mPG(Post-graduate),4:mGL(Graduate or lower)}
    \item{DES}{Desirable education level required for the principal job:
				(1:dPD(Postdoc),2:dARQ(Advanced Research Qualification),
				3:dPG(Post-graduate),4:dGL(Graduate or lower)}
    \item{PJREL}{Is your principal job related to your advanced research
				qualification degree:(1:H(High),2:M(Medium),3:L(Low))}
    \item{FOSAT}{Field of science and technology         
			(1:NS(Natural Sciences),2:ET(Engineering and technology),
			3:MH(Medical and health sciences),4:AS(Agricultural sciences),
			5:SS(Social Sciences),6:H(Humanities))}
    \item{SOF}{Principal source of financial support during your research studies
	          (1:F(Fellowship),2:T(Teaching),3:OE (Other Employment),4:R(Reimbursement)
			  ,5:LPSO (LoanPersonalSavingsOther)}	
  }
}
\source{
Spanish Statistical Office (Survey on Human Resources in Science and Technology, 2006):
http://www.ine.es/prodyser/micro_recurciencia.htm
}
\examples{
data(PhD_nomCyL)
}
\keyword{datasets}
