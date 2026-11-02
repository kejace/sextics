import Reflect.Core0
namespace HodgeComputer.Reflect
set_option maxRecDepth 100000
-- scaleT only: 100 scalings of a 100-term poly, summed lengths
def scaleAll (p q : SP) : Nat := p.foldr (fun t acc => (scaleT t q).length + acc) 0
example : scaleAll (mkPoly 1 100) (mkPoly 2 100) = 10000 := by decide +kernel
end HodgeComputer.Reflect
