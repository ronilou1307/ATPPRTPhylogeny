{-# LANGUAGE ExtendedDefaultRules #-}
module Main where
import Bio.Alignment
import Bio.Alphabet
import Bio.Sequence
import IModel
import MCMC
import Probability
import SModel
import Tree
import Tree.Newick
import qualified Data.IntMap as IntMap
import qualified Data.JSON as J
import qualified Data.Text.IO as T
import Probability.Logger
import System.Environment
import System.FilePath

sample_smodel  alpha = do {pi <- sample (symmetric_dirichlet_on (letters alpha) 1)
;alpha_2 <- sample (log_laplace 6 2)
;let {result = ((SModel.lg alpha +> SModel.plus_f' alpha pi) +> unit_mixture) +> SModel.gamma_rates alpha_2 4}
;let {loggers = ["f:pi" %=% pi, "Rates.gamma:alpha" %=% alpha_2]}
;return (result, loggers)
}

sample_imodel  topology = do {rate <- sample (log_laplace (negate 4) 0.707)
;mean_length <- sample (shifted_exponential 10 1)
;let {result = IModel.rs07 rate mean_length topology}
;let {loggers = ["rs07:rate" %=% rate, "rs07:mean_length" %=% mean_length]}
;return (result, loggers)
}

sample_scale  = sample (gamma 0.5 2)

sample_branch_lengths tree = sample (iidMap (getUEdgesSet tree) (gamma 0.5 (2 / intToDouble (numBranches tree))))

sample_topology taxa = uniform_labelled_topology taxa

model sequenceData logParamsTSV logTree [logA] = do {let {taxa = getTaxa sequenceData}
;
;topology <- RanSamplingRate 0 (sample_topology taxa)
;branch_lengths <- RanSamplingRate 0 (sample_branch_lengths topology)
;let {tree = branch_length_tree topology branch_lengths}
;RanSamplingRate 1 (PerformTKEffect (add_tree_moves tree))
;let {tlength = tree_length tree}
;scale1 <- sample_scale
;RanSamplingRate 2 (PerformTKEffect (add_move (scale_means_only_slice [scale1] (IntMap.elems branch_lengths))))
;RanSamplingRate 1 (PerformTKEffect (add_move (scale_means_only_MH [scale1] (IntMap.elems branch_lengths))))
;(smodel, log_smodel) <- sample_smodel aa
;(imodel, log_imodel) <- sample_imodel tree
;
;let {sequence_lengths = get_sequence_lengths sequenceData}
;(alignment, properties_A) <- sampleWithProps (phyloAlignment tree imodel scale1 sequence_lengths)
;properties <- observe sequenceData (phyloCTMC tree alignment smodel scale1)
;
;let {alignment_length = alignmentLength alignment}
;let {num_indels = totalNumIndels alignment}
;let {total_length_indels = totalLengthIndels alignment}
;let {prior_A = ln (probability properties_A)}
;let {anc_alignment = toFasta (prop_anc_seqs properties)}
;let {substs = prop_n_muts properties}
;let {p1_loggers = ["|A|" %=% alignment_length, "#indels" %=% num_indels, "|indels|" %=% total_length_indels, "prior_A" %=% prior_A, "likelihood" %=% ln (prop_likelihood properties), "#substs" %=% substs]}
;
;let {alignmentLengths = [alignment_length]}
;let {scale = scale1}
;let {loggers = ["|T|" %=% tlength, "scale1" %=% scale1, "scale1*|T|" %=% (scale1 * tlength), "S1" %>% log_smodel, "I1" %>% log_imodel, "P1" %>% p1_loggers, "scale" %=% scale, "scale*|T|" %=% (scale * tlength), "|A|" %=% sum alignmentLengths, "#indels" %=% sum [num_indels], "|indels|" %=% sum [total_length_indels], "#substs" %=% sum [substs], "prior_A" %=% sum [prior_A]]}
;
;
;addLogger $ logParamsTSV loggers
;
;addLogger $ logTree (addInternalLabels (scale_branch_lengths scale tree))
;
;addLogger $ (every 10 $ logA anc_alignment)
;
;return loggers
}

main = do {[directory] <- getArgs
;sequenceData <- mkUnalignedCharacterData aa <$> load_sequences "/home/rlouden/phylogeny/HisGS_21May.txt"
;
;logParamsTSV <- tsvLogger (directory </> "C1.log") ["iter"]
;
;logTree <- treeLogger (directory </> "C1.trees")
;
;logA <- alignmentLogger (directory </> "C1.P1.fastas")
;
;mymodel <- makeMCMCModel $ model sequenceData logParamsTSV logTree [logA]
;
;runMCMC 200000 mymodel
}
