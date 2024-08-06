library(tidyverse)
library(foreach)
library(openai)
library(ontologyIndex)
library(ontologySimilarity)

# Need this to access OpenAI API
# If you don't have one, you need to register an account with OpenAI
Sys.setenv(OPENAI_API_KEY = "put your key here")


# load ontology graph
TO    <- ontologyIndex::get_ontology(file = "ontology/to.obo")
TO.ic <- ontologySimilarity::descendants_IC(TO)

