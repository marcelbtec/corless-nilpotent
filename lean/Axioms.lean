import CorlessNilpotent

-- Corless's Conjecture 5.1, in his own terms (vectors `u` with entries in `{1, i}`)
#print axioms Corless.corless_conjecture
-- the number of nilpotent matrices in dimension `2^(k+1) - 1`
#print axioms Corless.card_nilpotent
-- the same classification for an arbitrary `u : ℕ → ℂ` and for the matching polynomial
#print axioms Corless.T_nilpotent_iff
#print axioms Corless.corless
