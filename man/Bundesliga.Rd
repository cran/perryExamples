\name{Bundesliga}
\alias{Bundesliga}
\docType{data}
\title{
Austrian Bundesliga football player data
}
\description{
The data set contains information on the market value of midfielders and
forwards in the Austrian \emph{Bundesliga}, together with player information
and performance measures.  The data are collected for the (still ongoing)
2013/14 season, with performance measures referring to competitions on the
Austrian level (Bundesliga, Cup) or the European level (UEFA Champions League,
UEFA Europa League).  Only players with complete information are included in
the data set.
}
\usage{data("Bundesliga")}
\format{
  A data frame with 123 observations on the following 20 variables.
  \describe{
    \item{\code{Player}}{factor; the player's name.}
    \item{\code{Team}}{factor; the player's team.}
    \item{\code{MarketValue}}{numeric; the player's market value (in Euros).}
    \item{\code{Age}}{numeric; the player's age (in years).}
    \item{\code{Height}}{numeric; the player's height (in cm).}
    \item{\code{Foreign}}{a dummy variable indicating whether the player is
    foreign or Austrian.}
    \item{\code{Forward}}{a dummy variable indicating whether the player is a
    forward or midfielder.}
    \item{\code{BothFeet}}{a dummy variable indicating whether the player is
     equally strong with both feet or has one stronger foot.}
    \item{\code{AtClub}}{numeric; the number of seasons the player is with his
    current club (at the upcoming transfer window).}
    \item{\code{Contract}}{numeric; the remaining number of seasons in the
    player's contract (at the upcoming transfer window).}
    \item{\code{Matches}}{numeric; the number of matches in which the player
    was on the field.}
    \item{\code{Goals}}{numeric; the number of goals the player scored.}
    \item{\code{OwnGoals}}{numeric; the number of own goals the player scored.}
    \item{\code{Assists}}{numeric; the number of assists the player gave.}
    \item{\code{Yellow}}{numeric; the number of yellow cards the player
    received.}
    \item{\code{YellowRed}}{numeric; the number of times the player was sent
    off with two yellow cards within one game.}
    \item{\code{Red}}{numeric; the number of times the player was sent off with
    a red card.}
    \item{\code{SubOn}}{numeric; the number of times the player was substituted
    on.}
    \item{\code{SubOff}}{numeric; the number of times the player was
    substituted off.}
    \item{\code{Minutes}}{numeric; the total number of minutes the player was
    on the field.}
  }
}
\source{
The data were scraped from \code{http://www.transfermarkt.com} on 2014-03-02.
}
\examples{
data("Bundesliga")
summary(Bundesliga)
}
\keyword{datasets}
