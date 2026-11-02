import Reflect.Kron

set_option maxHeartbeats 4000000
set_option maxRecDepth 8000

set_option linter.all false
set_option profiler true
set_option profiler.threshold 1000

namespace HodgeComputer.Instances.KronWitness

open MvPolynomial HodgeTestbed HodgeComputer.Kron

noncomputable def tF : PT := .node (.node (.node (.node (.leaf ⟨0, 0, 0, 0, 0, 6⟩ (-2)) (.leaf ⟨6, 0, 0, 0, 0, 0⟩ (-116))) (.node (.leaf ⟨4, 2, 0, 0, 0, 0⟩ (132)) (.node (.leaf ⟨2, 4, 0, 0, 0, 0⟩ (14)) (.leaf ⟨0, 6, 0, 0, 0, 0⟩ (-30))))) (.node (.node (.leaf ⟨0, 0, 2, 0, 0, 4⟩ (40)) (.node (.leaf ⟨2, 0, 2, 0, 0, 2⟩ (1)) (.leaf ⟨0, 0, 4, 0, 0, 2⟩ (80)))) (.node (.leaf ⟨0, 0, 6, 0, 0, 0⟩ (22)) (.node (.leaf ⟨4, 0, 0, 2, 0, 0⟩ (-171)) (.leaf ⟨2, 2, 0, 2, 0, 0⟩ (-52)))))) (.node (.node (.node (.leaf ⟨0, 4, 0, 2, 0, 0⟩ (62)) (.leaf ⟨2, 0, 0, 4, 0, 0⟩ (33))) (.node (.leaf ⟨0, 2, 0, 4, 0, 0⟩ (-34)) (.node (.leaf ⟨0, 0, 0, 6, 0, 0⟩ (-3)) (.leaf ⟨0, 0, 0, 0, 2, 4⟩ (8))))) (.node (.node (.leaf ⟨0, 0, 2, 0, 2, 2⟩ (6)) (.node (.leaf ⟨0, 0, 4, 0, 2, 0⟩ (-56)) (.leaf ⟨0, 2, 0, 2, 2, 0⟩ (1)))) (.node (.leaf ⟨0, 0, 0, 0, 4, 2⟩ (40)) (.node (.leaf ⟨0, 0, 2, 0, 4, 0⟩ (32)) (.leaf ⟨0, 0, 0, 0, 6, 0⟩ (2))))))

noncomputable def tA : PT := .node (.leaf ⟨3, 1, 1, 0, 1, 0⟩ (4)) (.node (.leaf ⟨1, 3, 1, 0, 1, 0⟩ (5)) (.leaf ⟨1, 1, 1, 2, 1, 0⟩ (-5)))

noncomputable def tD : PT := .node (.node (.node (.node (.node (.node (.node (.leaf ⟨0, 0, 0, 0, 0, 10⟩ (-24)) (.leaf ⟨2, 0, 0, 0, 0, 8⟩ (-2782080))) (.node (.leaf ⟨4, 0, 0, 0, 0, 6⟩ (43776)) (.node (.leaf ⟨6, 0, 0, 0, 0, 4⟩ (360)) (.leaf ⟨8, 0, 0, 0, 0, 2⟩ (21434112))))) (.node (.node (.leaf ⟨0, 2, 0, 0, 0, 8⟩ (2782080)) (.node (.leaf ⟨2, 2, 0, 0, 0, 6⟩ (-87552)) (.leaf ⟨4, 2, 0, 0, 0, 4⟩ (-1080)))) (.node (.leaf ⟨6, 2, 0, 0, 0, 2⟩ (-85736448)) (.node (.leaf ⟨0, 4, 0, 0, 0, 6⟩ (43776)) (.leaf ⟨2, 4, 0, 0, 0, 4⟩ (1080)))))) (.node (.node (.node (.leaf ⟨4, 4, 0, 0, 0, 2⟩ (128604672)) (.node (.leaf ⟨0, 6, 0, 0, 0, 4⟩ (-360)) (.leaf ⟨2, 6, 0, 0, 0, 2⟩ (-85736448)))) (.node (.leaf ⟨0, 8, 0, 0, 0, 2⟩ (21434112)) (.node (.leaf ⟨0, 0, 2, 0, 0, 8⟩ (-96)) (.leaf ⟨2, 0, 2, 0, 0, 6⟩ (3642240))))) (.node (.node (.leaf ⟨4, 0, 2, 0, 0, 4⟩ (-78432)) (.node (.leaf ⟨8, 0, 2, 0, 0, 0⟩ (-2381568)) (.leaf ⟨0, 2, 2, 0, 0, 6⟩ (-3642240)))) (.node (.leaf ⟨2, 2, 2, 0, 0, 4⟩ (156864)) (.node (.leaf ⟨6, 2, 2, 0, 0, 0⟩ (9526272)) (.leaf ⟨0, 4, 2, 0, 0, 4⟩ (-78432))))))) (.node (.node (.node (.node (.leaf ⟨4, 4, 2, 0, 0, 0⟩ (-14289408)) (.node (.leaf ⟨2, 6, 2, 0, 0, 0⟩ (9526272)) (.leaf ⟨0, 8, 2, 0, 0, 0⟩ (-2381568)))) (.node (.leaf ⟨0, 0, 4, 0, 0, 6⟩ (480)) (.node (.leaf ⟨2, 0, 4, 0, 0, 4⟩ (8033760)) (.leaf ⟨4, 0, 4, 0, 0, 2⟩ (34656))))) (.node (.node (.leaf ⟨0, 2, 4, 0, 0, 4⟩ (-8033760)) (.node (.leaf ⟨2, 2, 4, 0, 0, 2⟩ (-69312)) (.leaf ⟨0, 4, 4, 0, 0, 2⟩ (34656)))) (.node (.leaf ⟨0, 0, 6, 0, 0, 4⟩ (-24)) (.node (.leaf ⟨2, 0, 6, 0, 0, 2⟩ (2761920)) (.leaf ⟨0, 2, 6, 0, 0, 2⟩ (-2761920)))))) (.node (.node (.node (.leaf ⟨2, 0, 8, 0, 0, 0⟩ (-1118880)) (.node (.leaf ⟨0, 2, 8, 0, 0, 0⟩ (1118880)) (.leaf ⟨0, 0, 0, 2, 0, 8⟩ (-6398832)))) (.node (.leaf ⟨2, 0, 0, 2, 0, 6⟩ (41472)) (.node (.leaf ⟨4, 0, 0, 2, 0, 4⟩ (125160)) (.leaf ⟨6, 0, 0, 2, 0, 2⟩ (72007296))))) (.node (.node (.leaf ⟨0, 2, 0, 2, 0, 6⟩ (-41472)) (.node (.leaf ⟨2, 2, 0, 2, 0, 4⟩ (-250320)) (.leaf ⟨4, 2, 0, 2, 0, 2⟩ (-216021888)))) (.node (.leaf ⟨0, 4, 0, 2, 0, 4⟩ (125160)) (.node (.leaf ⟨2, 4, 0, 2, 0, 2⟩ (216021888)) (.leaf ⟨0, 6, 0, 2, 0, 2⟩ (-72007296)))))))) (.node (.node (.node (.node (.node (.leaf ⟨0, 0, 2, 2, 0, 6⟩ (8377536)) (.node (.leaf ⟨2, 0, 2, 2, 0, 4⟩ (-74292)) (.leaf ⟨4, 0, 2, 2, 0, 2⟩ (430272)))) (.node (.leaf ⟨6, 0, 2, 2, 0, 0⟩ (-8000832)) (.node (.leaf ⟨0, 2, 2, 2, 0, 4⟩ (74292)) (.leaf ⟨2, 2, 2, 2, 0, 2⟩ (-860544))))) (.node (.node (.leaf ⟨4, 2, 2, 2, 0, 0⟩ (24002496)) (.node (.leaf ⟨0, 4, 2, 2, 0, 2⟩ (430272)) (.leaf ⟨2, 4, 2, 2, 0, 0⟩ (-24002496)))) (.node (.leaf ⟨0, 6, 2, 2, 0, 0⟩ (8000832)) (.node (.leaf ⟨0, 0, 4, 2, 0, 4⟩ (18477576)) (.leaf ⟨2, 0, 4, 2, 0, 2⟩ (32832)))))) (.node (.node (.node (.leaf ⟨4, 0, 4, 2, 0, 0⟩ (25920)) (.node (.leaf ⟨0, 2, 4, 2, 0, 2⟩ (-32832)) (.leaf ⟨2, 2, 4, 2, 0, 0⟩ (-51840)))) (.node (.leaf ⟨0, 4, 4, 2, 0, 0⟩ (25920)) (.node (.leaf ⟨0, 0, 6, 2, 0, 2⟩ (6352824)) (.leaf ⟨0, 0, 8, 2, 0, 0⟩ (-2573424))))) (.node (.node (.leaf ⟨0, 0, 0, 4, 0, 6⟩ (552)) (.node (.leaf ⟨2, 0, 0, 4, 0, 4⟩ (227352)) (.leaf ⟨4, 0, 0, 4, 0, 2⟩ (82834032)))) (.node (.leaf ⟨6, 0, 0, 4, 0, 0⟩ (-1392)) (.node (.leaf ⟨0, 2, 0, 4, 0, 4⟩ (-227352)) (.leaf ⟨2, 2, 0, 4, 0, 2⟩ (-165668064))))))) (.node (.node (.node (.node (.leaf ⟨4, 2, 0, 4, 0, 0⟩ (4176)) (.node (.leaf ⟨0, 4, 0, 4, 0, 2⟩ (82834032)) (.leaf ⟨2, 4, 0, 4, 0, 0⟩ (-4176)))) (.node (.leaf ⟨0, 6, 0, 4, 0, 0⟩ (1392)) (.node (.leaf ⟨0, 0, 2, 4, 0, 4⟩ (-552)) (.leaf ⟨2, 0, 2, 4, 0, 2⟩ (784860))))) (.node (.node (.leaf ⟨4, 0, 2, 4, 0, 0⟩ (-9203712)) (.node (.leaf ⟨0, 2, 2, 4, 0, 2⟩ (-784860)) (.leaf ⟨2, 2, 2, 4, 0, 0⟩ (18407424)))) (.node (.leaf ⟨0, 4, 2, 4, 0, 0⟩ (-9203712)) (.node (.leaf ⟨0, 0, 4, 4, 0, 2⟩ (1416)) (.leaf ⟨2, 0, 4, 4, 0, 0⟩ (47280)))))) (.node (.node (.node (.leaf ⟨0, 2, 4, 4, 0, 0⟩ (-47280)) (.node (.leaf ⟨0, 0, 6, 4, 0, 0⟩ (264)) (.leaf ⟨0, 0, 0, 6, 0, 4⟩ (540)))) (.node (.leaf ⟨2, 0, 0, 6, 0, 2⟩ (32055204)) (.node (.leaf ⟨4, 0, 0, 6, 0, 0⟩ (-2052)) (.leaf ⟨0, 2, 0, 6, 0, 2⟩ (-32055204))))) (.node (.node (.leaf ⟨2, 2, 0, 6, 0, 0⟩ (4104)) (.node (.leaf ⟨0, 4, 0, 6, 0, 0⟩ (-2052)) (.leaf ⟨0, 0, 2, 6, 0, 2⟩ (1992)))) (.node (.leaf ⟨2, 0, 2, 6, 0, 0⟩ (-3561600)) (.node (.leaf ⟨0, 2, 2, 6, 0, 0⟩ (3561600)) (.leaf ⟨0, 0, 4, 6, 0, 0⟩ (120))))))))) (.node (.node (.node (.node (.node (.node (.leaf ⟨0, 0, 0, 8, 0, 2⟩ (-2919744)) (.leaf ⟨2, 0, 0, 8, 0, 0⟩ (396))) (.node (.leaf ⟨0, 2, 0, 8, 0, 0⟩ (-396)) (.node (.leaf ⟨0, 0, 2, 8, 0, 0⟩ (324408)) (.leaf ⟨0, 0, 0, 10, 0, 0⟩ (-36))))) (.node (.node (.leaf ⟨0, 0, 0, 0, 2, 8⟩ (96)) (.node (.leaf ⟨2, 0, 0, 0, 2, 6⟩ (-3642240)) (.leaf ⟨4, 0, 0, 0, 2, 4⟩ (78432)))) (.node (.leaf ⟨8, 0, 0, 0, 2, 0⟩ (2381568)) (.node (.leaf ⟨0, 2, 0, 0, 2, 6⟩ (3642240)) (.leaf ⟨2, 2, 0, 0, 2, 4⟩ (-156864)))))) (.node (.node (.node (.leaf ⟨6, 2, 0, 0, 2, 0⟩ (-9526272)) (.node (.leaf ⟨0, 4, 0, 0, 2, 4⟩ (78432)) (.leaf ⟨4, 4, 0, 0, 2, 0⟩ (14289408)))) (.node (.leaf ⟨2, 6, 0, 0, 2, 0⟩ (-9526272)) (.node (.leaf ⟨0, 8, 0, 0, 2, 0⟩ (2381568)) (.leaf ⟨0, 0, 2, 0, 2, 6⟩ (-960))))) (.node (.node (.leaf ⟨2, 0, 2, 0, 2, 4⟩ (-16067520)) (.node (.leaf ⟨4, 0, 2, 0, 2, 2⟩ (-69312)) (.leaf ⟨0, 2, 2, 0, 2, 4⟩ (16067520)))) (.node (.leaf ⟨2, 2, 2, 0, 2, 2⟩ (138624)) (.node (.leaf ⟨0, 4, 2, 0, 2, 2⟩ (-69312)) (.leaf ⟨0, 0, 4, 0, 2, 4⟩ (72))))))) (.node (.node (.node (.node (.leaf ⟨2, 0, 4, 0, 2, 2⟩ (-8285760)) (.node (.leaf ⟨0, 2, 4, 0, 2, 2⟩ (8285760)) (.leaf ⟨2, 0, 6, 0, 2, 0⟩ (4475520)))) (.node (.leaf ⟨0, 2, 6, 0, 2, 0⟩ (-4475520)) (.node (.leaf ⟨0, 0, 0, 2, 2, 6⟩ (-8377536)) (.leaf ⟨2, 0, 0, 2, 2, 4⟩ (74292))))) (.node (.node (.leaf ⟨4, 0, 0, 2, 2, 2⟩ (-430272)) (.node (.leaf ⟨6, 0, 0, 2, 2, 0⟩ (8000832)) (.leaf ⟨0, 2, 0, 2, 2, 4⟩ (-74292)))) (.node (.leaf ⟨2, 2, 0, 2, 2, 2⟩ (860544)) (.node (.leaf ⟨4, 2, 0, 2, 2, 0⟩ (-24002496)) (.leaf ⟨0, 4, 0, 2, 2, 2⟩ (-430272)))))) (.node (.node (.node (.leaf ⟨2, 4, 0, 2, 2, 0⟩ (24002496)) (.node (.leaf ⟨0, 6, 0, 2, 2, 0⟩ (-8000832)) (.leaf ⟨0, 0, 2, 2, 2, 4⟩ (-36955152)))) (.node (.leaf ⟨2, 0, 2, 2, 2, 2⟩ (-65664)) (.node (.leaf ⟨4, 0, 2, 2, 2, 0⟩ (-51840)) (.leaf ⟨0, 2, 2, 2, 2, 2⟩ (65664))))) (.node (.node (.leaf ⟨2, 2, 2, 2, 2, 0⟩ (103680)) (.node (.leaf ⟨0, 4, 2, 2, 2, 0⟩ (-51840)) (.leaf ⟨0, 0, 4, 2, 2, 2⟩ (-19058472)))) (.node (.leaf ⟨0, 0, 6, 2, 2, 0⟩ (10293696)) (.node (.leaf ⟨0, 0, 0, 4, 2, 4⟩ (552)) (.leaf ⟨2, 0, 0, 4, 2, 2⟩ (-784860)))))))) (.node (.node (.node (.node (.node (.leaf ⟨4, 0, 0, 4, 2, 0⟩ (9203712)) (.node (.leaf ⟨0, 2, 0, 4, 2, 2⟩ (784860)) (.leaf ⟨2, 2, 0, 4, 2, 0⟩ (-18407424)))) (.node (.leaf ⟨0, 4, 0, 4, 2, 0⟩ (9203712)) (.node (.leaf ⟨0, 0, 2, 4, 2, 2⟩ (-2832)) (.leaf ⟨2, 0, 2, 4, 2, 0⟩ (-94560))))) (.node (.node (.leaf ⟨0, 2, 2, 4, 2, 0⟩ (94560)) (.node (.leaf ⟨0, 0, 4, 4, 2, 0⟩ (-792)) (.leaf ⟨0, 0, 0, 6, 2, 2⟩ (-1992)))) (.node (.leaf ⟨2, 0, 0, 6, 2, 0⟩ (3561600)) (.node (.leaf ⟨0, 2, 0, 6, 2, 0⟩ (-3561600)) (.leaf ⟨0, 0, 2, 6, 2, 0⟩ (-240)))))) (.node (.node (.node (.leaf ⟨0, 0, 0, 8, 2, 0⟩ (-324408)) (.node (.leaf ⟨0, 0, 0, 0, 4, 6⟩ (480)) (.leaf ⟨2, 0, 0, 0, 4, 4⟩ (8033760)))) (.node (.leaf ⟨4, 0, 0, 0, 4, 2⟩ (34656)) (.node (.leaf ⟨0, 2, 0, 0, 4, 4⟩ (-8033760)) (.leaf ⟨2, 2, 0, 0, 4, 2⟩ (-69312))))) (.node (.node (.leaf ⟨0, 4, 0, 0, 4, 2⟩ (34656)) (.node (.leaf ⟨0, 0, 2, 0, 4, 4⟩ (-72)) (.leaf ⟨2, 0, 2, 0, 4, 2⟩ (8285760)))) (.node (.leaf ⟨0, 2, 2, 0, 4, 2⟩ (-8285760)) (.node (.leaf ⟨2, 0, 4, 0, 4, 0⟩ (-6713280)) (.leaf ⟨0, 2, 4, 0, 4, 0⟩ (6713280))))))) (.node (.node (.node (.node (.leaf ⟨0, 0, 0, 2, 4, 4⟩ (18477576)) (.node (.leaf ⟨2, 0, 0, 2, 4, 2⟩ (32832)) (.leaf ⟨4, 0, 0, 2, 4, 0⟩ (25920)))) (.node (.leaf ⟨0, 2, 0, 2, 4, 2⟩ (-32832)) (.node (.leaf ⟨2, 2, 0, 2, 4, 0⟩ (-51840)) (.leaf ⟨0, 4, 0, 2, 4, 0⟩ (25920))))) (.node (.node (.leaf ⟨0, 0, 2, 2, 4, 2⟩ (19058472)) (.node (.leaf ⟨0, 0, 4, 2, 4, 0⟩ (-15440544)) (.leaf ⟨0, 0, 0, 4, 4, 2⟩ (1416)))) (.node (.leaf ⟨2, 0, 0, 4, 4, 0⟩ (47280)) (.node (.leaf ⟨0, 2, 0, 4, 4, 0⟩ (-47280)) (.leaf ⟨0, 0, 2, 4, 4, 0⟩ (792)))))) (.node (.node (.node (.leaf ⟨0, 0, 0, 6, 4, 0⟩ (120)) (.node (.leaf ⟨0, 0, 0, 0, 6, 4⟩ (24)) (.leaf ⟨2, 0, 0, 0, 6, 2⟩ (-2761920)))) (.node (.leaf ⟨0, 2, 0, 0, 6, 2⟩ (2761920)) (.node (.leaf ⟨2, 0, 2, 0, 6, 0⟩ (4475520)) (.leaf ⟨0, 2, 2, 0, 6, 0⟩ (-4475520))))) (.node (.node (.leaf ⟨0, 0, 0, 2, 6, 2⟩ (-6352824)) (.node (.leaf ⟨0, 0, 2, 2, 6, 0⟩ (10293696)) (.leaf ⟨0, 0, 0, 4, 6, 0⟩ (-264)))) (.node (.leaf ⟨2, 0, 0, 0, 8, 0⟩ (-1118880)) (.node (.leaf ⟨0, 2, 0, 0, 8, 0⟩ (1118880)) (.leaf ⟨0, 0, 0, 2, 8, 0⟩ (-2573424)))))))))

noncomputable def tQ : PT := .leaf ⟨14, 0, 0, 0, 0, 0⟩ 0

noncomputable def tN0 : PT := .node (.node (.node (.node (.node (.node (.leaf ⟨0, 1, 1, 0, 1, 8⟩ (89424)) (.leaf ⟨2, 1, 1, 0, 1, 6⟩ (1152))) (.node (.leaf ⟨4, 1, 1, 0, 1, 4⟩ (5)) (.node (.leaf ⟨6, 1, 1, 0, 1, 2⟩ (227808)) (.leaf ⟨0, 3, 1, 0, 1, 6⟩ (-720))))) (.node (.node (.leaf ⟨2, 3, 1, 0, 1, 4⟩ (-5)) (.node (.leaf ⟨4, 3, 1, 0, 1, 2⟩ (-334656)) (.leaf ⟨2, 5, 1, 0, 1, 2⟩ (-14112)))) (.node (.leaf ⟨0, 7, 1, 0, 1, 2⟩ (120960)) (.node (.leaf ⟨0, 1, 3, 0, 1, 6⟩ (-117072)) (.leaf ⟨2, 1, 3, 0, 1, 4⟩ (-2064)))))) (.node (.node (.node (.leaf ⟨6, 1, 3, 0, 1, 0⟩ (-25312)) (.leaf ⟨0, 3, 3, 0, 1, 4⟩ (1290))) (.node (.leaf ⟨4, 3, 3, 0, 1, 0⟩ (37184)) (.node (.leaf ⟨2, 5, 3, 0, 1, 0⟩ (1568)) (.leaf ⟨0, 7, 3, 0, 1, 0⟩ (-13440))))) (.node (.node (.leaf ⟨0, 1, 5, 0, 1, 4⟩ (-258228)) (.node (.leaf ⟨2, 1, 5, 0, 1, 2⟩ (912)) (.leaf ⟨0, 3, 5, 0, 1, 2⟩ (-570)))) (.node (.leaf ⟨0, 1, 7, 0, 1, 2⟩ (-88776)) (.node (.leaf ⟨0, 1, 9, 0, 1, 0⟩ (35964)) (.leaf ⟨0, 1, 1, 2, 1, 6⟩ (720))))))) (.node (.node (.node (.node (.leaf ⟨2, 1, 1, 2, 1, 4⟩ (-1435)) (.leaf ⟨4, 1, 1, 2, 1, 2⟩ (461146))) (.node (.leaf ⟨0, 3, 1, 2, 1, 4⟩ (1872)) (.node (.leaf ⟨2, 3, 1, 2, 1, 2⟩ (-176875)) (.leaf ⟨0, 5, 1, 2, 1, 2⟩ (-284271))))) (.node (.node (.leaf ⟨0, 1, 3, 2, 1, 4⟩ (-1290)) (.node (.leaf ⟨2, 1, 3, 2, 1, 2⟩ (-4980)) (.leaf ⟨4, 1, 3, 2, 1, 0⟩ (-51240)))) (.node (.leaf ⟨0, 3, 3, 2, 1, 2⟩ (6474)) (.node (.leaf ⟨2, 3, 3, 2, 1, 0⟩ (19656)) (.leaf ⟨0, 5, 3, 2, 1, 0⟩ (31584)))))) (.node (.node (.node (.leaf ⟨0, 1, 5, 2, 1, 2⟩ (570)) (.leaf ⟨2, 1, 5, 2, 1, 0⟩ (-300))) (.node (.leaf ⟨0, 3, 5, 2, 1, 0⟩ (390)) (.node (.leaf ⟨0, 1, 1, 4, 1, 4⟩ (720)) (.leaf ⟨2, 1, 1, 4, 1, 2⟩ (221251))))) (.node (.node (.leaf ⟨4, 1, 1, 4, 1, 0⟩ (8)) (.node (.leaf ⟨0, 3, 1, 4, 1, 2⟩ (186747)) (.leaf ⟨2, 3, 1, 4, 1, 0⟩ (-20)))) (.node (.leaf ⟨0, 5, 1, 4, 1, 0⟩ (12)) (.node (.leaf ⟨0, 1, 3, 4, 1, 2⟩ (2490)) (.leaf ⟨2, 1, 3, 4, 1, 0⟩ (-24584)))))))) (.node (.node (.node (.node (.node (.leaf ⟨0, 3, 3, 4, 1, 0⟩ (-20748)) (.leaf ⟨0, 1, 5, 4, 1, 0⟩ (150))) (.node (.leaf ⟨0, 1, 1, 6, 1, 2⟩ (44604)) (.node (.leaf ⟨2, 1, 1, 6, 1, 0⟩ (-10)) (.leaf ⟨0, 3, 1, 6, 1, 0⟩ (15))))) (.node (.node (.leaf ⟨0, 1, 3, 6, 1, 0⟩ (-4956)) (.node (.leaf ⟨0, 1, 1, 0, 3, 6⟩ (117072)) (.leaf ⟨2, 1, 1, 0, 3, 4⟩ (2064)))) (.node (.leaf ⟨6, 1, 1, 0, 3, 0⟩ (25312)) (.node (.leaf ⟨0, 3, 1, 0, 3, 4⟩ (-1290)) (.leaf ⟨4, 3, 1, 0, 3, 0⟩ (-37184)))))) (.node (.node (.node (.leaf ⟨2, 5, 1, 0, 3, 0⟩ (-1568)) (.leaf ⟨0, 7, 1, 0, 3, 0⟩ (13440))) (.node (.leaf ⟨0, 1, 3, 0, 3, 4⟩ (516456)) (.node (.leaf ⟨2, 1, 3, 0, 3, 2⟩ (-1824)) (.leaf ⟨0, 3, 3, 0, 3, 2⟩ (1140))))) (.node (.node (.leaf ⟨0, 1, 5, 0, 3, 2⟩ (266328)) (.node (.leaf ⟨0, 1, 7, 0, 3, 0⟩ (-143856)) (.leaf ⟨0, 1, 1, 2, 3, 4⟩ (1290)))) (.node (.leaf ⟨2, 1, 1, 2, 3, 2⟩ (4980)) (.node (.leaf ⟨4, 1, 1, 2, 3, 0⟩ (51240)) (.leaf ⟨0, 3, 1, 2, 3, 2⟩ (-6474))))))) (.node (.node (.node (.node (.leaf ⟨2, 3, 1, 2, 3, 0⟩ (-19656)) (.leaf ⟨0, 5, 1, 2, 3, 0⟩ (-31584))) (.node (.leaf ⟨0, 1, 3, 2, 3, 2⟩ (-1140)) (.node (.leaf ⟨2, 1, 3, 2, 3, 0⟩ (600)) (.leaf ⟨0, 3, 3, 2, 3, 0⟩ (-780))))) (.node (.node (.leaf ⟨0, 1, 1, 4, 3, 2⟩ (-2490)) (.node (.leaf ⟨2, 1, 1, 4, 3, 0⟩ (24584)) (.leaf ⟨0, 3, 1, 4, 3, 0⟩ (20748)))) (.node (.leaf ⟨0, 1, 3, 4, 3, 0⟩ (-300)) (.node (.leaf ⟨0, 1, 1, 6, 3, 0⟩ (4956)) (.leaf ⟨0, 1, 1, 0, 5, 4⟩ (-258228)))))) (.node (.node (.node (.leaf ⟨2, 1, 1, 0, 5, 2⟩ (912)) (.leaf ⟨0, 3, 1, 0, 5, 2⟩ (-570))) (.node (.leaf ⟨0, 1, 3, 0, 5, 2⟩ (-266328)) (.node (.leaf ⟨0, 1, 5, 0, 5, 0⟩ (215784)) (.leaf ⟨0, 1, 1, 2, 5, 2⟩ (570))))) (.node (.node (.leaf ⟨2, 1, 1, 2, 5, 0⟩ (-300)) (.node (.leaf ⟨0, 3, 1, 2, 5, 0⟩ (390)) (.leaf ⟨0, 1, 1, 4, 5, 0⟩ (150)))) (.node (.leaf ⟨0, 1, 1, 0, 7, 2⟩ (88776)) (.node (.leaf ⟨0, 1, 3, 0, 7, 0⟩ (-143856)) (.leaf ⟨0, 1, 1, 0, 9, 0⟩ (35964))))))))

noncomputable def tN1 : PT := .node (.node (.node (.node (.node (.node (.leaf ⟨1, 0, 1, 0, 1, 8⟩ (89424)) (.leaf ⟨3, 0, 1, 0, 1, 6⟩ (2016))) (.node (.leaf ⟨5, 0, 1, 0, 1, 4⟩ (15)) (.node (.leaf ⟨7, 0, 1, 0, 1, 2⟩ (925344)) (.leaf ⟨1, 2, 1, 0, 1, 6⟩ (-1584))))) (.node (.node (.leaf ⟨3, 2, 1, 0, 1, 4⟩ (-25)) (.leaf ⟨5, 2, 1, 0, 1, 2⟩ (-2427264))) (.node (.leaf ⟨1, 4, 1, 0, 1, 4⟩ (10)) (.node (.leaf ⟨3, 4, 1, 0, 1, 2⟩ (2078496)) (.leaf ⟨1, 6, 1, 0, 1, 2⟩ (-576576)))))) (.node (.node (.node (.leaf ⟨1, 0, 3, 0, 1, 6⟩ (-117072)) (.leaf ⟨3, 0, 3, 0, 1, 4⟩ (-3612))) (.node (.leaf ⟨7, 0, 3, 0, 1, 0⟩ (-102816)) (.node (.leaf ⟨1, 2, 3, 0, 1, 4⟩ (2838)) (.leaf ⟨5, 2, 3, 0, 1, 0⟩ (269696))))) (.node (.node (.leaf ⟨3, 4, 3, 0, 1, 0⟩ (-230944)) (.node (.leaf ⟨1, 6, 3, 0, 1, 0⟩ (64064)) (.leaf ⟨1, 0, 5, 0, 1, 4⟩ (-258228)))) (.node (.leaf ⟨3, 0, 5, 0, 1, 2⟩ (1596)) (.node (.leaf ⟨1, 2, 5, 0, 1, 2⟩ (-1254)) (.leaf ⟨1, 0, 7, 0, 1, 2⟩ (-88776))))))) (.node (.node (.node (.node (.leaf ⟨1, 0, 9, 0, 1, 0⟩ (35964)) (.leaf ⟨1, 0, 1, 2, 1, 6⟩ (720))) (.node (.leaf ⟨3, 0, 1, 2, 1, 4⟩ (-561)) (.node (.leaf ⟨5, 0, 1, 2, 1, 2⟩ (1951980)) (.leaf ⟨1, 2, 1, 2, 1, 4⟩ (998))))) (.node (.node (.leaf ⟨3, 2, 1, 2, 1, 2⟩ (-3158543)) (.node (.leaf ⟨1, 4, 1, 2, 1, 2⟩ (1206563)) (.leaf ⟨1, 0, 3, 2, 1, 4⟩ (-1290)))) (.node (.leaf ⟨3, 0, 3, 2, 1, 2⟩ (-1992)) (.node (.leaf ⟨5, 0, 3, 2, 1, 0⟩ (-216888)) (.leaf ⟨1, 2, 3, 2, 1, 2⟩ (3486)))))) (.node (.node (.node (.leaf ⟨3, 2, 3, 2, 1, 0⟩ (350952)) (.leaf ⟨1, 4, 3, 2, 1, 0⟩ (-134064))) (.node (.leaf ⟨1, 0, 5, 2, 1, 2⟩ (570)) (.node (.leaf ⟨3, 0, 5, 2, 1, 0⟩ (-120)) (.leaf ⟨1, 2, 5, 2, 1, 0⟩ (210))))) (.node (.node (.leaf ⟨1, 0, 1, 4, 1, 4⟩ (720)) (.node (.leaf ⟨3, 0, 1, 4, 1, 2⟩ (1037247)) (.leaf ⟨1, 2, 1, 4, 1, 2⟩ (-629249)))) (.node (.leaf ⟨3, 2, 1, 4, 1, 0⟩ (-4)) (.node (.leaf ⟨1, 4, 1, 4, 1, 0⟩ (4)) (.leaf ⟨1, 0, 3, 4, 1, 2⟩ (2490)))))))) (.node (.node (.node (.node (.node (.leaf ⟨3, 0, 3, 4, 1, 0⟩ (-115248)) (.leaf ⟨1, 2, 3, 4, 1, 0⟩ (69916))) (.node (.leaf ⟨1, 0, 5, 4, 1, 0⟩ (150)) (.node (.leaf ⟨1, 0, 1, 6, 1, 2⟩ (44604)) (.leaf ⟨1, 2, 1, 6, 1, 0⟩ (5))))) (.node (.node (.leaf ⟨1, 0, 3, 6, 1, 0⟩ (-4956)) (.node (.leaf ⟨1, 0, 1, 0, 3, 6⟩ (117072)) (.leaf ⟨3, 0, 1, 0, 3, 4⟩ (3612)))) (.node (.leaf ⟨7, 0, 1, 0, 3, 0⟩ (102816)) (.node (.leaf ⟨1, 2, 1, 0, 3, 4⟩ (-2838)) (.leaf ⟨5, 2, 1, 0, 3, 0⟩ (-269696)))))) (.node (.node (.node (.leaf ⟨3, 4, 1, 0, 3, 0⟩ (230944)) (.leaf ⟨1, 6, 1, 0, 3, 0⟩ (-64064))) (.node (.leaf ⟨1, 0, 3, 0, 3, 4⟩ (516456)) (.node (.leaf ⟨3, 0, 3, 0, 3, 2⟩ (-3192)) (.leaf ⟨1, 2, 3, 0, 3, 2⟩ (2508))))) (.node (.node (.leaf ⟨1, 0, 5, 0, 3, 2⟩ (266328)) (.node (.leaf ⟨1, 0, 7, 0, 3, 0⟩ (-143856)) (.leaf ⟨1, 0, 1, 2, 3, 4⟩ (1290)))) (.node (.leaf ⟨3, 0, 1, 2, 3, 2⟩ (1992)) (.node (.leaf ⟨5, 0, 1, 2, 3, 0⟩ (216888)) (.leaf ⟨1, 2, 1, 2, 3, 2⟩ (-3486))))))) (.node (.node (.node (.node (.leaf ⟨3, 2, 1, 2, 3, 0⟩ (-350952)) (.leaf ⟨1, 4, 1, 2, 3, 0⟩ (134064))) (.node (.leaf ⟨1, 0, 3, 2, 3, 2⟩ (-1140)) (.node (.leaf ⟨3, 0, 3, 2, 3, 0⟩ (240)) (.leaf ⟨1, 2, 3, 2, 3, 0⟩ (-420))))) (.node (.node (.leaf ⟨1, 0, 1, 4, 3, 2⟩ (-2490)) (.node (.leaf ⟨3, 0, 1, 4, 3, 0⟩ (115248)) (.leaf ⟨1, 2, 1, 4, 3, 0⟩ (-69916)))) (.node (.leaf ⟨1, 0, 3, 4, 3, 0⟩ (-300)) (.node (.leaf ⟨1, 0, 1, 6, 3, 0⟩ (4956)) (.leaf ⟨1, 0, 1, 0, 5, 4⟩ (-258228)))))) (.node (.node (.node (.leaf ⟨3, 0, 1, 0, 5, 2⟩ (1596)) (.leaf ⟨1, 2, 1, 0, 5, 2⟩ (-1254))) (.node (.leaf ⟨1, 0, 3, 0, 5, 2⟩ (-266328)) (.node (.leaf ⟨1, 0, 5, 0, 5, 0⟩ (215784)) (.leaf ⟨1, 0, 1, 2, 5, 2⟩ (570))))) (.node (.node (.leaf ⟨3, 0, 1, 2, 5, 0⟩ (-120)) (.node (.leaf ⟨1, 2, 1, 2, 5, 0⟩ (210)) (.leaf ⟨1, 0, 1, 4, 5, 0⟩ (150)))) (.node (.leaf ⟨1, 0, 1, 0, 7, 2⟩ (88776)) (.node (.leaf ⟨1, 0, 3, 0, 7, 0⟩ (-143856)) (.leaf ⟨1, 0, 1, 0, 9, 0⟩ (35964))))))))

noncomputable def tN2 : PT := .node (.node (.node (.node (.node (.node (.leaf ⟨1, 1, 0, 0, 1, 8⟩ (-108)) (.leaf ⟨5, 1, 0, 0, 1, 4⟩ (199296))) (.node (.leaf ⟨7, 1, 0, 0, 1, 2⟩ (480)) (.node (.leaf ⟨3, 3, 0, 0, 1, 4⟩ (-398592)) (.leaf ⟨5, 3, 0, 0, 1, 2⟩ (-1440))))) (.node (.node (.leaf ⟨1, 5, 0, 0, 1, 4⟩ (199296)) (.leaf ⟨3, 5, 0, 0, 1, 2⟩ (1440))) (.node (.leaf ⟨1, 7, 0, 0, 1, 2⟩ (-480)) (.node (.leaf ⟨3, 1, 2, 0, 1, 4⟩ (5)) (.leaf ⟨5, 1, 2, 0, 1, 2⟩ (-8304)))))) (.node (.node (.node (.leaf ⟨1, 3, 2, 0, 1, 4⟩ (-5)) (.leaf ⟨3, 3, 2, 0, 1, 2⟩ (16608))) (.node (.leaf ⟨1, 5, 2, 0, 1, 2⟩ (-8304)) (.node (.leaf ⟨1, 1, 4, 0, 1, 4⟩ (1386)) (.leaf ⟨5, 1, 4, 0, 1, 0⟩ (119024))))) (.node (.node (.leaf ⟨3, 3, 4, 0, 1, 0⟩ (-238048)) (.leaf ⟨1, 5, 4, 0, 1, 0⟩ (119024))) (.node (.leaf ⟨1, 1, 6, 0, 1, 2⟩ (234)) (.node (.leaf ⟨1, 1, 0, 2, 1, 6⟩ (-108)) (.leaf ⟨3, 1, 0, 2, 1, 4⟩ (425952))))))) (.node (.node (.node (.node (.leaf ⟨5, 1, 0, 2, 1, 2⟩ (1128)) (.leaf ⟨7, 1, 0, 2, 1, 0⟩ (-3672))) (.node (.leaf ⟨1, 3, 0, 2, 1, 4⟩ (-425952)) (.node (.leaf ⟨3, 3, 0, 2, 1, 2⟩ (-2256)) (.leaf ⟨5, 3, 0, 2, 1, 0⟩ (11016))))) (.node (.node (.leaf ⟨1, 5, 0, 2, 1, 2⟩ (1128)) (.leaf ⟨3, 5, 0, 2, 1, 0⟩ (-11016))) (.node (.leaf ⟨1, 7, 0, 2, 1, 0⟩ (3672)) (.node (.leaf ⟨1, 1, 2, 2, 1, 4⟩ (1301)) (.leaf ⟨3, 1, 2, 2, 1, 2⟩ (-17708)))))) (.node (.node (.node (.leaf ⟨1, 3, 2, 2, 1, 2⟩ (17708)) (.leaf ⟨1, 1, 4, 2, 1, 2⟩ (1332))) (.node (.leaf ⟨3, 1, 4, 2, 1, 0⟩ (254388)) (.node (.leaf ⟨1, 3, 4, 2, 1, 0⟩ (-254388)) (.leaf ⟨1, 1, 6, 2, 1, 0⟩ (1008))))) (.node (.node (.leaf ⟨1, 1, 0, 4, 1, 4⟩ (233136)) (.node (.leaf ⟨3, 1, 0, 4, 1, 2⟩ (741)) (.leaf ⟨5, 1, 0, 4, 1, 0⟩ (-7746)))) (.node (.leaf ⟨1, 3, 0, 4, 1, 2⟩ (-741)) (.node (.leaf ⟨3, 3, 0, 4, 1, 0⟩ (15492)) (.leaf ⟨1, 5, 0, 4, 1, 0⟩ (-7746)))))))) (.node (.node (.node (.node (.node (.leaf ⟨1, 1, 2, 4, 1, 2⟩ (-9719)) (.leaf ⟨3, 1, 2, 4, 1, 0⟩ (8))) (.node (.leaf ⟨1, 3, 2, 4, 1, 0⟩ (-8)) (.node (.leaf ⟨1, 1, 4, 4, 1, 0⟩ (139234)) (.leaf ⟨1, 1, 0, 6, 1, 2⟩ (-177))))) (.node (.node (.leaf ⟨3, 1, 0, 6, 1, 0⟩ (-4116)) (.leaf ⟨1, 3, 0, 6, 1, 0⟩ (4116))) (.node (.leaf ⟨1, 1, 2, 6, 1, 0⟩ (-10)) (.node (.leaf ⟨1, 1, 0, 8, 1, 0⟩ (-177)) (.leaf ⟨1, 1, 0, 0, 3, 6⟩ (432)))))) (.node (.node (.node (.leaf ⟨5, 1, 0, 0, 3, 2⟩ (357072)) (.leaf ⟨3, 3, 0, 0, 3, 2⟩ (-714144))) (.node (.leaf ⟨1, 5, 0, 0, 3, 2⟩ (357072)) (.node (.leaf ⟨1, 1, 2, 0, 3, 4⟩ (-3546)) (.leaf ⟨5, 1, 2, 0, 3, 0⟩ (-276800))))) (.node (.node (.leaf ⟨3, 3, 2, 0, 3, 0⟩ (553600)) (.node (.leaf ⟨1, 5, 2, 0, 3, 0⟩ (-276800)) (.leaf ⟨1, 1, 4, 0, 3, 2⟩ (-360)))) (.node (.leaf ⟨1, 1, 0, 2, 3, 4⟩ (-864)) (.node (.leaf ⟨3, 1, 0, 2, 3, 2⟩ (763125)) (.leaf ⟨1, 3, 0, 2, 3, 2⟩ (-763125))))))) (.node (.node (.node (.node (.leaf ⟨1, 1, 2, 2, 3, 2⟩ (-1170)) (.leaf ⟨3, 1, 2, 2, 3, 0⟩ (-591600))) (.node (.leaf ⟨1, 3, 2, 2, 3, 0⟩ (591600)) (.node (.leaf ⟨1, 1, 4, 2, 3, 0⟩ (-2934)) (.leaf ⟨1, 1, 0, 4, 3, 2⟩ (417717))))) (.node (.node (.leaf ⟨3, 1, 0, 4, 3, 0⟩ (-12)) (.leaf ⟨1, 3, 0, 4, 3, 0⟩ (12))) (.node (.leaf ⟨1, 1, 2, 4, 3, 0⟩ (-323800)) (.node (.leaf ⟨1, 1, 0, 6, 3, 0⟩ (15)) (.leaf ⟨1, 1, 0, 0, 5, 4⟩ (2160)))))) (.node (.node (.node (.leaf ⟨5, 1, 0, 0, 5, 0⟩ (157776)) (.leaf ⟨3, 3, 0, 0, 5, 0⟩ (-315552))) (.node (.leaf ⟨1, 5, 0, 0, 5, 0⟩ (157776)) (.node (.leaf ⟨1, 1, 2, 0, 5, 2⟩ (18)) (.leaf ⟨1, 1, 0, 2, 5, 2⟩ (-162))))) (.node (.node (.leaf ⟨3, 1, 0, 2, 5, 0⟩ (337212)) (.node (.leaf ⟨1, 3, 0, 2, 5, 0⟩ (-337212)) (.leaf ⟨1, 1, 2, 2, 5, 0⟩ (2844)))) (.node (.leaf ⟨1, 1, 0, 4, 5, 0⟩ (184566)) (.node (.leaf ⟨1, 1, 0, 0, 7, 2⟩ (108)) (.leaf ⟨1, 1, 0, 2, 7, 0⟩ (-918))))))))

noncomputable def tN3 : PT := .node (.node (.node (.node (.node (.leaf ⟨1, 1, 1, 1, 1, 6⟩ (-864)) (.leaf ⟨3, 1, 1, 1, 1, 4⟩ (-10))) (.node (.leaf ⟨5, 1, 1, 1, 1, 2⟩ (-697536)) (.leaf ⟨1, 3, 1, 1, 1, 4⟩ (10)))) (.node (.node (.leaf ⟨3, 3, 1, 1, 1, 2⟩ (1395072)) (.leaf ⟨1, 5, 1, 1, 1, 2⟩ (-697536))) (.node (.leaf ⟨1, 1, 3, 1, 1, 4⟩ (1548)) (.leaf ⟨5, 1, 3, 1, 1, 0⟩ (77504))))) (.node (.node (.node (.leaf ⟨3, 3, 3, 1, 1, 0⟩ (-155008)) (.leaf ⟨1, 5, 3, 1, 1, 0⟩ (77504))) (.node (.leaf ⟨1, 1, 5, 1, 1, 2⟩ (-684)) (.leaf ⟨1, 1, 1, 3, 1, 4⟩ (-874)))) (.node (.node (.leaf ⟨3, 1, 1, 3, 1, 2⟩ (-1490834)) (.leaf ⟨1, 3, 1, 3, 1, 2⟩ (1490834))) (.node (.leaf ⟨1, 1, 3, 3, 1, 2⟩ (-2988)) (.node (.leaf ⟨3, 1, 3, 3, 1, 0⟩ (165648)) (.leaf ⟨1, 3, 3, 3, 1, 0⟩ (-165648))))))) (.node (.node (.node (.node (.leaf ⟨1, 1, 5, 3, 1, 0⟩ (-180)) (.leaf ⟨1, 1, 1, 5, 1, 2⟩ (-815996))) (.node (.leaf ⟨3, 1, 1, 5, 1, 0⟩ (8)) (.leaf ⟨1, 3, 1, 5, 1, 0⟩ (-8)))) (.node (.node (.leaf ⟨1, 1, 3, 5, 1, 0⟩ (90664)) (.leaf ⟨1, 1, 1, 7, 1, 0⟩ (-10))) (.node (.leaf ⟨1, 1, 1, 1, 3, 4⟩ (-1548)) (.node (.leaf ⟨5, 1, 1, 1, 3, 0⟩ (-77504)) (.leaf ⟨3, 3, 1, 1, 3, 0⟩ (155008)))))) (.node (.node (.node (.leaf ⟨1, 5, 1, 1, 3, 0⟩ (-77504)) (.leaf ⟨1, 1, 3, 1, 3, 2⟩ (1368))) (.node (.leaf ⟨1, 1, 1, 3, 3, 2⟩ (2988)) (.leaf ⟨3, 1, 1, 3, 3, 0⟩ (-165648)))) (.node (.node (.leaf ⟨1, 3, 1, 3, 3, 0⟩ (165648)) (.leaf ⟨1, 1, 3, 3, 3, 0⟩ (360))) (.node (.leaf ⟨1, 1, 1, 5, 3, 0⟩ (-90664)) (.node (.leaf ⟨1, 1, 1, 1, 5, 2⟩ (-684)) (.leaf ⟨1, 1, 1, 3, 5, 0⟩ (-180)))))))

noncomputable def tN4 : PT := .node (.node (.node (.node (.node (.node (.leaf ⟨1, 1, 1, 0, 0, 8⟩ (-108)) (.leaf ⟨5, 1, 1, 0, 0, 4⟩ (199296))) (.node (.leaf ⟨7, 1, 1, 0, 0, 2⟩ (480)) (.node (.leaf ⟨3, 3, 1, 0, 0, 4⟩ (-398592)) (.leaf ⟨5, 3, 1, 0, 0, 2⟩ (-1440))))) (.node (.node (.leaf ⟨1, 5, 1, 0, 0, 4⟩ (199296)) (.leaf ⟨3, 5, 1, 0, 0, 2⟩ (1440))) (.node (.leaf ⟨1, 7, 1, 0, 0, 2⟩ (-480)) (.node (.leaf ⟨1, 1, 3, 0, 0, 6⟩ (864)) (.leaf ⟨3, 1, 3, 0, 0, 4⟩ (15)))))) (.node (.node (.node (.leaf ⟨5, 1, 3, 0, 0, 2⟩ (689232)) (.leaf ⟨1, 3, 3, 0, 0, 4⟩ (-15))) (.node (.leaf ⟨3, 3, 3, 0, 0, 2⟩ (-1378464)) (.node (.leaf ⟨1, 5, 3, 0, 0, 2⟩ (689232)) (.leaf ⟨1, 1, 5, 0, 0, 4⟩ (-162))))) (.node (.node (.leaf ⟨5, 1, 5, 0, 0, 0⟩ (41520)) (.leaf ⟨3, 3, 5, 0, 0, 0⟩ (-83040))) (.node (.leaf ⟨1, 5, 5, 0, 0, 0⟩ (41520)) (.node (.leaf ⟨1, 1, 7, 0, 0, 2⟩ (918)) (.leaf ⟨1, 1, 1, 2, 0, 6⟩ (-108))))))) (.node (.node (.node (.node (.leaf ⟨3, 1, 1, 2, 0, 4⟩ (425952)) (.leaf ⟨5, 1, 1, 2, 0, 2⟩ (1128))) (.node (.leaf ⟨7, 1, 1, 2, 0, 0⟩ (-3672)) (.node (.leaf ⟨1, 3, 1, 2, 0, 4⟩ (-425952)) (.leaf ⟨3, 3, 1, 2, 0, 2⟩ (-2256))))) (.node (.node (.leaf ⟨5, 3, 1, 2, 0, 0⟩ (11016)) (.leaf ⟨1, 5, 1, 2, 0, 2⟩ (1128))) (.node (.leaf ⟨3, 5, 1, 2, 0, 0⟩ (-11016)) (.node (.leaf ⟨1, 7, 1, 2, 0, 0⟩ (3672)) (.leaf ⟨1, 1, 3, 2, 0, 4⟩ (2175)))))) (.node (.node (.node (.leaf ⟨3, 1, 3, 2, 0, 2⟩ (1473126)) (.leaf ⟨1, 3, 3, 2, 0, 2⟩ (-1473126))) (.node (.leaf ⟨1, 1, 5, 2, 0, 2⟩ (4320)) (.node (.leaf ⟨3, 1, 5, 2, 0, 0⟩ (88740)) (.leaf ⟨1, 3, 5, 2, 0, 0⟩ (-88740))))) (.node (.node (.leaf ⟨1, 1, 7, 2, 0, 0⟩ (1188)) (.node (.leaf ⟨1, 1, 1, 4, 0, 4⟩ (233136)) (.leaf ⟨3, 1, 1, 4, 0, 2⟩ (741)))) (.node (.leaf ⟨5, 1, 1, 4, 0, 0⟩ (-7746)) (.node (.leaf ⟨1, 3, 1, 4, 0, 2⟩ (-741)) (.leaf ⟨3, 3, 1, 4, 0, 0⟩ (15492)))))))) (.node (.node (.node (.node (.node (.leaf ⟨1, 5, 1, 4, 0, 0⟩ (-7746)) (.leaf ⟨1, 1, 3, 4, 0, 2⟩ (806277))) (.node (.leaf ⟨1, 1, 5, 4, 0, 0⟩ (48570)) (.node (.leaf ⟨1, 1, 1, 6, 0, 2⟩ (-177)) (.leaf ⟨3, 1, 1, 6, 0, 0⟩ (-4116))))) (.node (.node (.leaf ⟨1, 3, 1, 6, 0, 0⟩ (4116)) (.leaf ⟨1, 1, 1, 8, 0, 0⟩ (-177))) (.node (.leaf ⟨1, 1, 1, 0, 2, 6⟩ (-432)) (.node (.leaf ⟨3, 1, 1, 0, 2, 4⟩ (-10)) (.leaf ⟨5, 1, 1, 0, 2, 2⟩ (-340464)))))) (.node (.node (.node (.leaf ⟨1, 3, 1, 0, 2, 4⟩ (10)) (.leaf ⟨3, 3, 1, 0, 2, 2⟩ (680928))) (.node (.leaf ⟨1, 5, 1, 0, 2, 2⟩ (-340464)) (.node (.leaf ⟨1, 1, 3, 0, 2, 4⟩ (-450)) (.leaf ⟨5, 1, 3, 0, 2, 0⟩ (-121792))))) (.node (.node (.leaf ⟨3, 3, 3, 0, 2, 0⟩ (243584)) (.node (.leaf ⟨1, 5, 3, 0, 2, 0⟩ (-121792)) (.leaf ⟨1, 1, 5, 0, 2, 2⟩ (-2412)))) (.node (.leaf ⟨1, 1, 1, 2, 2, 4⟩ (-1738)) (.node (.leaf ⟨3, 1, 1, 2, 2, 2⟩ (-727709)) (.leaf ⟨1, 3, 1, 2, 2, 2⟩ (727709))))))) (.node (.node (.node (.node (.leaf ⟨1, 1, 3, 2, 2, 2⟩ (-7146)) (.leaf ⟨3, 1, 3, 2, 2, 0⟩ (-260304))) (.node (.leaf ⟨1, 3, 3, 2, 2, 0⟩ (260304)) (.node (.leaf ⟨1, 1, 5, 2, 2, 0⟩ (-3474)) (.leaf ⟨1, 1, 1, 4, 2, 2⟩ (-398279))))) (.node (.node (.leaf ⟨3, 1, 1, 4, 2, 0⟩ (-4)) (.leaf ⟨1, 3, 1, 4, 2, 0⟩ (4))) (.node (.leaf ⟨1, 1, 3, 4, 2, 0⟩ (-142472)) (.node (.leaf ⟨1, 1, 1, 6, 2, 0⟩ (5)) (.leaf ⟨1, 1, 1, 0, 4, 4⟩ (612)))))) (.node (.node (.node (.leaf ⟨5, 1, 1, 0, 4, 0⟩ (80272)) (.leaf ⟨3, 3, 1, 0, 4, 0⟩ (-160544))) (.node (.leaf ⟨1, 5, 1, 0, 4, 0⟩ (80272)) (.node (.leaf ⟨1, 1, 3, 0, 4, 2⟩ (2070)) (.leaf ⟨1, 1, 1, 2, 4, 2⟩ (2826))))) (.node (.node (.leaf ⟨3, 1, 1, 2, 4, 0⟩ (171564)) (.node (.leaf ⟨1, 3, 1, 2, 4, 0⟩ (-171564)) (.leaf ⟨1, 1, 3, 2, 4, 0⟩ (3384)))) (.node (.leaf ⟨1, 1, 1, 4, 4, 0⟩ (93902)) (.node (.leaf ⟨1, 1, 1, 0, 6, 2⟩ (-576)) (.leaf ⟨1, 1, 1, 2, 6, 0⟩ (-1098))))))))

noncomputable def tN5 : PT := .node (.node (.node (.node (.node (.leaf ⟨1, 1, 1, 0, 1, 7⟩ (-864)) (.leaf ⟨3, 1, 1, 0, 1, 5⟩ (-10))) (.node (.leaf ⟨5, 1, 1, 0, 1, 3⟩ (-697536)) (.leaf ⟨1, 3, 1, 0, 1, 5⟩ (10)))) (.node (.node (.leaf ⟨3, 3, 1, 0, 1, 3⟩ (1395072)) (.leaf ⟨1, 5, 1, 0, 1, 3⟩ (-697536))) (.node (.leaf ⟨1, 1, 3, 0, 1, 5⟩ (1548)) (.leaf ⟨5, 1, 3, 0, 1, 1⟩ (77504))))) (.node (.node (.node (.leaf ⟨3, 3, 3, 0, 1, 1⟩ (-155008)) (.leaf ⟨1, 5, 3, 0, 1, 1⟩ (77504))) (.node (.leaf ⟨1, 1, 5, 0, 1, 3⟩ (-684)) (.leaf ⟨1, 1, 1, 2, 1, 5⟩ (-874)))) (.node (.node (.leaf ⟨3, 1, 1, 2, 1, 3⟩ (-1490834)) (.leaf ⟨1, 3, 1, 2, 1, 3⟩ (1490834))) (.node (.leaf ⟨1, 1, 3, 2, 1, 3⟩ (-2988)) (.node (.leaf ⟨3, 1, 3, 2, 1, 1⟩ (165648)) (.leaf ⟨1, 3, 3, 2, 1, 1⟩ (-165648))))))) (.node (.node (.node (.node (.leaf ⟨1, 1, 5, 2, 1, 1⟩ (-180)) (.leaf ⟨1, 1, 1, 4, 1, 3⟩ (-815996))) (.node (.leaf ⟨3, 1, 1, 4, 1, 1⟩ (8)) (.leaf ⟨1, 3, 1, 4, 1, 1⟩ (-8)))) (.node (.node (.leaf ⟨1, 1, 3, 4, 1, 1⟩ (90664)) (.leaf ⟨1, 1, 1, 6, 1, 1⟩ (-10))) (.node (.leaf ⟨1, 1, 1, 0, 3, 5⟩ (-1548)) (.node (.leaf ⟨5, 1, 1, 0, 3, 1⟩ (-77504)) (.leaf ⟨3, 3, 1, 0, 3, 1⟩ (155008)))))) (.node (.node (.node (.leaf ⟨1, 5, 1, 0, 3, 1⟩ (-77504)) (.leaf ⟨1, 1, 3, 0, 3, 3⟩ (1368))) (.node (.leaf ⟨1, 1, 1, 2, 3, 3⟩ (2988)) (.leaf ⟨3, 1, 1, 2, 3, 1⟩ (-165648)))) (.node (.node (.leaf ⟨1, 3, 1, 2, 3, 1⟩ (165648)) (.leaf ⟨1, 1, 3, 2, 3, 1⟩ (360))) (.node (.leaf ⟨1, 1, 1, 4, 3, 1⟩ (-90664)) (.node (.leaf ⟨1, 1, 1, 0, 5, 3⟩ (-684)) (.leaf ⟨1, 1, 1, 2, 5, 1⟩ (-180)))))))

noncomputable section

def F : Poly := toC tF
def A : Poly := toC tA

def D : Poly := toC tD
def Q : Poly := toC tQ

def N0 : Poly := toC tN0

def N1 : Poly := toC tN1

def N2 : Poly := toC tN2

def N3 : Poly := toC tN3

def N4 : Poly := toC tN4

def N5 : Poly := toC tN5

def N : Fin 6 → Poly := fun i => if i.val = 0 then N0 else if i.val = 1 then N1 else if i.val = 2 then N2 else if i.val = 3 then N3 else if i.val = 4 then N4 else N5

noncomputable def L1 : E := (.smulN 1 (.add (.mul (.pd 0 tF) (.poly tN0)) (.add (.mul (.pd 1 tF) (.poly tN1)) (.add (.mul (.pd 2 tF) (.poly tN2)) (.add (.mul (.pd 3 tF) (.poly tN3)) (.add (.mul (.pd 4 tF) (.poly tN4)) (.mul (.pd 5 tF) (.poly tN5))))))))
noncomputable def R1 : E := (.mul (.poly tA) (.poly tD))

theorem e1_core : denC L1 = denC R1 := soundN (m := 21) (W := 61) L1 R1 (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel)

theorem e1 : (1 : Poly) * directional F N = A*D := by
  simp only [directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]
  change (1:Poly) * (pderiv (0 : Fin 6) F * N0 + (pderiv (1 : Fin 6) F * N1 + (pderiv (2 : Fin 6) F * N2 + (pderiv (3 : Fin 6) F * N3 + (pderiv (4 : Fin 6) F * N4 + (pderiv (5 : Fin 6) F * N5)))))) = A*D
  have h := e1_core
  simp only [L1, R1, denC_smulN, denC_add, denC_mul, denC_poly, denC_pd, Nat.cast_one, Nat.cast_ofNat] at h
  simpa only [F, A, D, N0, N1, N2, N3, N4, N5] using h

noncomputable def L2 : E := (.sub (.mul (.poly tD) (.add (.pd 0 tN0) (.add (.pd 1 tN1) (.add (.pd 2 tN2) (.add (.pd 3 tN3) (.add (.pd 4 tN4) (.pd 5 tN5))))))) (.add (.mul (.pd 0 tD) (.poly tN0)) (.add (.mul (.pd 1 tD) (.poly tN1)) (.add (.mul (.pd 2 tD) (.poly tN2)) (.add (.mul (.pd 3 tD) (.poly tN3)) (.add (.mul (.pd 4 tD) (.poly tN4)) (.mul (.pd 5 tD) (.poly tN5))))))))
noncomputable def R2 : E := (.mul (.pow (.poly tF) 1) (.poly tQ))

theorem e2_core : denC L2 = denC R2 := soundN (m := 21) (W := 61) L2 R2 (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel)

theorem e2 : divergenceNumerator D N = F^1*Q := by
  simp only [divergenceNumerator, divergence, directional, Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero]
  change D * (pderiv (0 : Fin 6) N0 + (pderiv (1 : Fin 6) N1 + (pderiv (2 : Fin 6) N2 + (pderiv (3 : Fin 6) N3 + (pderiv (4 : Fin 6) N4 + (pderiv (5 : Fin 6) N5)))))) - (pderiv (0 : Fin 6) D * N0 + (pderiv (1 : Fin 6) D * N1 + (pderiv (2 : Fin 6) D * N2 + (pderiv (3 : Fin 6) D * N3 + (pderiv (4 : Fin 6) D * N4 + (pderiv (5 : Fin 6) D * N5)))))) = F^1*Q
  have h := e2_core
  simp only [L2, R2, denC_sub, denC_add, denC_mul, denC_pow, denC_poly, denC_pd] at h
  simpa only [F, D, Q, N0, N1, N2, N3, N4, N5] using h

theorem denominator_homogeneous : D.IsHomogeneous 10 := toC_isHomogeneous 10 tD (by decide +kernel)

theorem N0_homogeneous : N0.IsHomogeneous 11 := toC_isHomogeneous 11 tN0 (by decide +kernel)

theorem N1_homogeneous : N1.IsHomogeneous 11 := toC_isHomogeneous 11 tN1 (by decide +kernel)

theorem N2_homogeneous : N2.IsHomogeneous 11 := toC_isHomogeneous 11 tN2 (by decide +kernel)

theorem N3_homogeneous : N3.IsHomogeneous 11 := toC_isHomogeneous 11 tN3 (by decide +kernel)

theorem N4_homogeneous : N4.IsHomogeneous 11 := toC_isHomogeneous 11 tN4 (by decide +kernel)

theorem N5_homogeneous : N5.IsHomogeneous 11 := toC_isHomogeneous 11 tN5 (by decide +kernel)

theorem numerator_homogeneous : ∀ i, (N i).IsHomogeneous 11 := by
  intro i; fin_cases i
  · exact N0_homogeneous
  · exact N1_homogeneous
  · exact N2_homogeneous
  · exact N3_homogeneous
  · exact N4_homogeneous
  · exact N5_homogeneous

theorem denominator_ne_zero : D ≠ 0 := toC_ne_zero_of_evalPT ![1,1,1,1,1,1] tD (by decide +kernel)

def toUni : Poly →+* Polynomial Coeff := MvPolynomial.eval₂Hom (Polynomial.C.comp MvPolynomial.C) (fun i => if i.val = 0 then (Polynomial.X) else if i.val = 1 then (Polynomial.C (MvPolynomial.X 0)) else if i.val = 2 then (Polynomial.C (MvPolynomial.X 1)) else if i.val = 3 then (Polynomial.C (MvPolynomial.X 2)) else if i.val = 4 then (Polynomial.C (MvPolynomial.X 3)) else (Polynomial.C (MvPolynomial.X 4)))

def fromUni : Polynomial Coeff →+* Poly := Polynomial.eval₂RingHom (MvPolynomial.eval₂Hom MvPolynomial.C (fun i => if i.val = 0 then (MvPolynomial.X 1) else if i.val = 1 then (MvPolynomial.X 2) else if i.val = 2 then (MvPolynomial.X 3) else if i.val = 3 then (MvPolynomial.X 4) else (MvPolynomial.X 5))) (MvPolynomial.X 0)

@[simp] theorem toUni_X0 : toUni (MvPolynomial.X 0) = (Polynomial.X) := by simp only [toUni,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem toUni_X1 : toUni (MvPolynomial.X 1) = (Polynomial.C (MvPolynomial.X 0)) := by simp only [toUni,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem toUni_X2 : toUni (MvPolynomial.X 2) = (Polynomial.C (MvPolynomial.X 1)) := by simp only [toUni,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem toUni_X3 : toUni (MvPolynomial.X 3) = (Polynomial.C (MvPolynomial.X 2)) := by simp only [toUni,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem toUni_X4 : toUni (MvPolynomial.X 4) = (Polynomial.C (MvPolynomial.X 3)) := by simp only [toUni,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem toUni_X5 : toUni (MvPolynomial.X 5) = (Polynomial.C (MvPolynomial.X 4)) := by simp only [toUni,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem toUni_C (c : ℂ) : toUni (MvPolynomial.C c) = Polynomial.C (MvPolynomial.C c) := by simp [toUni]

@[simp] theorem fromUni_w0 : fromUni (Polynomial.C (MvPolynomial.X 0)) = (MvPolynomial.X 1) := by simp only [fromUni,Polynomial.coe_eval₂RingHom,Polynomial.eval₂_C,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem fromUni_w1 : fromUni (Polynomial.C (MvPolynomial.X 1)) = (MvPolynomial.X 2) := by simp only [fromUni,Polynomial.coe_eval₂RingHom,Polynomial.eval₂_C,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem fromUni_w2 : fromUni (Polynomial.C (MvPolynomial.X 2)) = (MvPolynomial.X 3) := by simp only [fromUni,Polynomial.coe_eval₂RingHom,Polynomial.eval₂_C,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem fromUni_w3 : fromUni (Polynomial.C (MvPolynomial.X 3)) = (MvPolynomial.X 4) := by simp only [fromUni,Polynomial.coe_eval₂RingHom,Polynomial.eval₂_C,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem fromUni_w4 : fromUni (Polynomial.C (MvPolynomial.X 4)) = (MvPolynomial.X 5) := by simp only [fromUni,Polynomial.coe_eval₂RingHom,Polynomial.eval₂_C,MvPolynomial.eval₂Hom_X']; rfl

@[simp] theorem fromUni_t : fromUni Polynomial.X = MvPolynomial.X 0 := by simp [fromUni]

theorem retract : fromUni.comp toUni = RingHom.id Poly := by
  apply MvPolynomial.ringHom_ext
  · intro c; simp [toUni,fromUni]
  · intro i; fin_cases i <;> simp only [RingHom.comp_apply,RingHom.id_apply]
    · change fromUni (toUni (X 0)) = X 0; simp only [toUni_X0,map_add,fromUni_t,fromUni_w0,fromUni_w1,fromUni_w2,fromUni_w3,fromUni_w4] <;> ring
    · change fromUni (toUni (X 1)) = X 1; simp only [toUni_X1,map_add,fromUni_t,fromUni_w0,fromUni_w1,fromUni_w2,fromUni_w3,fromUni_w4] <;> ring
    · change fromUni (toUni (X 2)) = X 2; simp only [toUni_X2,map_add,fromUni_t,fromUni_w0,fromUni_w1,fromUni_w2,fromUni_w3,fromUni_w4] <;> ring
    · change fromUni (toUni (X 3)) = X 3; simp only [toUni_X3,map_add,fromUni_t,fromUni_w0,fromUni_w1,fromUni_w2,fromUni_w3,fromUni_w4] <;> ring
    · change fromUni (toUni (X 4)) = X 4; simp only [toUni_X4,map_add,fromUni_t,fromUni_w0,fromUni_w1,fromUni_w2,fromUni_w3,fromUni_w4] <;> ring
    · change fromUni (toUni (X 5)) = X 5; simp only [toUni_X5,map_add,fromUni_t,fromUni_w0,fromUni_w1,fromUni_w2,fromUni_w3,fromUni_w4] <;> ring

def fUni : Polynomial Coeff := Polynomial.C ((MvPolynomial.C (-116 : ℂ))) * Polynomial.X^6 +
    Polynomial.C ((MvPolynomial.C (132 : ℂ)) * (MvPolynomial.X 0)^2 +
    (MvPolynomial.C (-171 : ℂ)) * (MvPolynomial.X 2)^2) * Polynomial.X^4 +
    Polynomial.C ((MvPolynomial.C (14 : ℂ)) * (MvPolynomial.X 0)^4 +
    (MvPolynomial.C (-52 : ℂ)) * (MvPolynomial.X 0)^2 * (MvPolynomial.X 2)^2 +
    (MvPolynomial.C (1 : ℂ)) * (MvPolynomial.X 1)^2 * (MvPolynomial.X 4)^2 +
    (MvPolynomial.C (33 : ℂ)) * (MvPolynomial.X 2)^4) * Polynomial.X^2 +
    Polynomial.C ((MvPolynomial.C (-30 : ℂ)) * (MvPolynomial.X 0)^6 +
    (MvPolynomial.C (62 : ℂ)) * (MvPolynomial.X 0)^4 * (MvPolynomial.X 2)^2 +
    (MvPolynomial.C (-34 : ℂ)) * (MvPolynomial.X 0)^2 * (MvPolynomial.X 2)^4 +
    (MvPolynomial.C (1 : ℂ)) * (MvPolynomial.X 0)^2 * (MvPolynomial.X 2)^2 * (MvPolynomial.X 3)^2 +
    (MvPolynomial.C (22 : ℂ)) * (MvPolynomial.X 1)^6 +
    (MvPolynomial.C (-56 : ℂ)) * (MvPolynomial.X 1)^4 * (MvPolynomial.X 3)^2 +
    (MvPolynomial.C (80 : ℂ)) * (MvPolynomial.X 1)^4 * (MvPolynomial.X 4)^2 +
    (MvPolynomial.C (32 : ℂ)) * (MvPolynomial.X 1)^2 * (MvPolynomial.X 3)^4 +
    (MvPolynomial.C (6 : ℂ)) * (MvPolynomial.X 1)^2 * (MvPolynomial.X 3)^2 * (MvPolynomial.X 4)^2 +
    (MvPolynomial.C (40 : ℂ)) * (MvPolynomial.X 1)^2 * (MvPolynomial.X 4)^4 +
    (MvPolynomial.C (-3 : ℂ)) * (MvPolynomial.X 2)^6 +
    (MvPolynomial.C (2 : ℂ)) * (MvPolynomial.X 3)^6 +
    (MvPolynomial.C (40 : ℂ)) * (MvPolynomial.X 3)^4 * (MvPolynomial.X 4)^2 +
    (MvPolynomial.C (8 : ℂ)) * (MvPolynomial.X 3)^2 * (MvPolynomial.X 4)^4 +
    (MvPolynomial.C (-2 : ℂ)) * (MvPolynomial.X 4)^6) * Polynomial.X^0

theorem toUni_F : toUni F = fUni := by
  simp only [F, tF, toC, fUni, map_add, map_mul, map_pow, toUni_C, toUni_X0,toUni_X1,toUni_X2,toUni_X3,toUni_X4,toUni_X5]
  push_cast
  polynomial

theorem fUni_degree : fUni.natDegree = 6 := by
  unfold fUni; compute_degree!

theorem fUni_lc : fUni.leadingCoeff = MvPolynomial.C (-116 : ℂ) := by
  rw [Polynomial.leadingCoeff, fUni_degree]
  simp only [fUni, Polynomial.coeff_add, Polynomial.coeff_C_mul_X_pow]
  norm_num

def point : Fin 5 → ℂ := fun i => if i.val = 0 then (1 : ℂ) else if i.val = 1 then (1 : ℂ) else if i.val = 2 then (1 : ℂ) else if i.val = 3 then (1 : ℂ) else (1 : ℂ)

def specialize : Coeff →+* ℂ := MvPolynomial.eval point

def fFiber : Polynomial ℂ := (Polynomial.C (-116 : ℂ)) * Polynomial.X^6 +
    (Polynomial.C (-39 : ℂ)) * Polynomial.X^4 +
    (Polynomial.C (-4 : ℂ)) * Polynomial.X^2 +
    (Polynomial.C (168 : ℂ))

def dFiber : Polynomial ℂ := (Polynomial.C (21434112 : ℂ)) * Polynomial.X^8 +
    (Polynomial.C (-13730184 : ℂ)) * Polynomial.X^6 +
    (Polynomial.C (-4413204 : ℂ)) * Polynomial.X^4 +
    (Polynomial.C (-6177144 : ℂ)) * Polynomial.X^2 +
    (Polynomial.C (-6431124 : ℂ))

def bezoutF : Polynomial ℂ := (Polynomial.C (-3885133032494643946112 / 911849364795315540749229 : ℂ)) * Polynomial.X^6 +
    (Polynomial.C (-54022649204988025936100 / 19148836660701626355733809 : ℂ)) * Polynomial.X^4 +
    (Polynomial.C (-9001019672329272828604 / 2735548094385946622247687 : ℂ)) * Polynomial.X^2 +
    (Polynomial.C (226069960164215706197401 / 57446509982104879067201427 : ℂ))

def bezoutD : Polynomial ℂ := (Polynomial.C (-3973929808914527174 / 172339529946314637201604281 : ℂ)) * Polynomial.X^4 +
    (Polynomial.C (-13025941903448747303 / 344679059892629274403208562 : ℂ)) * Polynomial.X^2 +
    (Polynomial.C (-12107840977419586639 / 229786039928419516268805708 : ℂ))

theorem fiber_F : (toUni F).map specialize = fFiber := by
  rw [toUni_F]
  norm_num [fUni,specialize,point,fFiber] <;> polynomial

theorem fiber_D : (toUni D).map specialize = dFiber := by
  simp only [D,tD,toC,map_add,map_mul,map_pow,toUni_C,toUni_X0,toUni_X1,toUni_X2,toUni_X3,toUni_X4,toUni_X5]
  push_cast
  norm_num [specialize,point,dFiber] <;> polynomial

theorem fiber_coprime : IsCoprime fFiber dFiber := by
  refine ⟨bezoutF, bezoutD, ?_⟩
  unfold bezoutF bezoutD fFiber dFiber
  polynomial

theorem denominator_legal : IsRelPrime F D := by
  apply relPrime_of_retract toUni fromUni retract
  apply relPrime_of_specialization specialize
  · rw [toUni_F, fUni_lc]
    exact (isUnit_iff_ne_zero.mpr (by norm_num : (-116 : ℂ) ≠ 0)).map MvPolynomial.C
  · rw [fiber_F, fiber_D]
    exact fiber_coprime.isRelPrime

theorem certificate : Verify F A 1 10 D N Q where
  denominator_ne_zero := denominator_ne_zero
  denominator_homogeneous := denominator_homogeneous
  numerator_homogeneous := numerator_homogeneous
  denominator_legal := denominator_legal
  identity := normalized_identity (by simpa using e1) e2

theorem solved : HasCertificate F A 1 := ⟨10, D, N, Q, certificate⟩

end
end HodgeComputer.Instances.KronWitness
