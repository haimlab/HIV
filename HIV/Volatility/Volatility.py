import tkinter
from tkinter import filedialog
import os
from ComputeStandardDevs import getVal
from Bio import SeqIO
import math
from numpy import linalg as LA
import numpy as np

def vals(name, place, el):
    if name == 'G':
        if place == 0:
            if el == 'A':
                return -0.8433408256
            elif el == 'C':
                return -0.8433408256
            elif el == 'D':
                return -0.9181034093
            elif el == 'E':
                return -0.639319686
            elif el == 'F':
                return -0.8433408256
            elif el == 'G':
                return -0.8433408256
            elif el == 'H':
                return -1.386786226
            elif el == 'I':
                return -0.8433408256
            elif el == 'K':
                return -1.024623924
            elif el == 'L':
                return -1.156683389
            elif el == 'M':
                return -0.8433408256
            elif el == 'N':
                return -1.053210671
            elif el == 'P':
                return -0.8433408256
            elif el == 'Q':
                return -0.8433408256
            elif el == 'R':
                return -1.205705075
            elif el == 'S':
                return -0.172860397
            elif el == 'T':
                return -0.9906088454
            elif el == 'V':
                return -0.8433408256
            elif el == 'W':
                return -0.8433408256
            elif el == 'Y':
                return -0.8433408256
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -2.2
            else:
                while True:
                    print('Broken ' + el)
                return -100
        if place == 1:
            if el == 'A':
                return -0.4640493953
            elif el == 'C':
                return -0.4640493953
            elif el == 'D':
                return -1.403275213
            elif el == 'E':
                return -0.7335853809
            elif el == 'F':
                return -0.4640493953
            elif el == 'G':
                return -0.4640493953
            elif el == 'H':
                return -0.776090809
            elif el == 'I':
                return 0.2940775637
            elif el == 'K':
                return -0.776090809
            elif el == 'L':
                return -0.4640493953
            elif el == 'M':
                return -0.4640493953
            elif el == 'N':
                return -0.7839769129
            elif el == 'P':
                return -0.4640493953
            elif el == 'Q':
                return -0.4640493953
            elif el == 'R':
                return -0.776090809
            elif el == 'S':
                return -0.7167583816
            elif el == 'T':
                return -0.6495398504
            elif el == 'V':
                return -0.4640493953
            elif el == 'W':
                return -0.4640493953
            elif el == 'Y':
                return -0.4640493953
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -2.2
            else:
                while True:
                    print('Broken ' + el)
                return -100
        if place == 2:
            if el == 'A':
                return -0.786797717
            elif el == 'C':
                return -0.786797717
            elif el == 'D':
                return -0.07329156295
            elif el == 'E':
                return -0.1360604223
            elif el == 'F':
                return -0.04847787444
            elif el == 'G':
                return -1.46648923
            elif el == 'H':
                return -0.4798840043
            elif el == 'I':
                return -1.32674822
            elif el == 'K':
                return -0.2562146956
            elif el == 'L':
                return -0.786797717
            elif el == 'M':
                return -0.786797717
            elif el == 'N':
                return -0.3081689069
            elif el == 'P':
                return -0.786797717
            elif el == 'Q':
                return -0.786797717
            elif el == 'R':
                return -0.703553313
            elif el == 'S':
                return -0.05240518592
            elif el == 'T':
                return -2.137180921
            elif el == 'V':
                return -0.786797717
            elif el == 'W':
                return -0.786797717
            elif el == 'Y':
                return -0.3054755433
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -2.2
            else:
                while True:
                    print('Broken ' + el)
                return -100
        if place == 3:
            if el == 'A':
                return -0.8448633008
            elif el == 'C':
                return -0.8448633008
            elif el == 'D':
                return -1.368665877
            elif el == 'E':
                return -1.368665877
            elif el == 'F':
                return -0.8448633008
            elif el == 'G':
                return -0.8448633008
            elif el == 'H':
                return -1.299082897
            elif el == 'I':
                return -0.8448633008
            elif el == 'K':
                return -1.299082897
            elif el == 'L':
                return -0.8448633008
            elif el == 'M':
                return -0.8448633008
            elif el == 'N':
                return -1.110448535
            elif el == 'P':
                return -0.8448633008
            elif el == 'Q':
                return -0.8448633008
            elif el == 'R':
                return -1.299082897
            elif el == 'S':
                return -0.6382200956
            elif el == 'T':
                return -0.7859212715
            elif el == 'V':
                return -0.8448633008
            elif el == 'W':
                return -0.8448633008
            elif el == 'Y':
                return -0.8448633008
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -2.2
            else:
                while True:
                    print('Broken ' + el)
                return -100
        if place == 4:
            if el == 'A':
                return -0.2614702126
            elif el == 'C':
                return -0.2614702126
            elif el == 'D':
                return -0.1675110553
            elif el == 'E':
                return -0.1675110553
            elif el == 'F':
                return -0.2614702126
            elif el == 'G':
                return -0.2614702126
            elif el == 'H':
                return -0.1675110553
            elif el == 'I':
                return -0.2614702126
            elif el == 'K':
                return -0.1675110553
            elif el == 'L':
                return -0.2614702126
            elif el == 'M':
                return -0.2614702126
            elif el == 'N':
                return -0.1369302597
            elif el == 'P':
                return -0.2614702126
            elif el == 'Q':
                return -0.6350900713
            elif el == 'R':
                return -0.1675110553
            elif el == 'S':
                return 0.2658180943
            elif el == 'T':
                return -0.5396786137
            elif el == 'V':
                return -0.2614702126
            elif el == 'W':
                return -0.2614702126
            elif el == 'Y':
                return -0.2614702126
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -2.2
            else:
                while True:
                    print('Broken ' + el)
                return -100
    if name == 'F':
        if place == 0:
            if el == 'A':
                return -0.3723497881
            elif el == 'C':
                return -0.4250787875
            elif el == 'D':
                return 0
            elif el == 'E':
                return 0
            elif el == 'F':
                return -0.4250787875
            elif el == 'G':
                return -0.5331361326
            elif el == 'H':
                return -0.5169944736
            elif el == 'I':
                return -0.4250787875
            elif el == 'K':
                return -0.5169944736
            elif el == 'L':
                return -0.4250787875
            elif el == 'M':
                return -0.4250787875
            elif el == 'N':
                return -0.3883523771
            elif el == 'P':
                return -0.4250787875
            elif el == 'Q':
                return -0.3697504419
            elif el == 'R':
                return -0.5169944736
            elif el == 'S':
                return -0.7567166597
            elif el == 'T':
                return -0.01998809453
            elif el == 'V':
                return -0.4250787875
            elif el == 'W':
                return -0.4250787875
            elif el == 'Y':
                return -0.4250787875
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -1.9
            else:
                while True:
                    print('Broken ' + el)
                return -100
        if place == 1:
            if el == 'A':
                return -0.5077579086
            elif el == 'C':
                return -0.5077579086
            elif el == 'D':
                return -0.5077579086
            elif el == 'E':
                return -0.5077579086
            elif el == 'F':
                return -0.5077579086
            elif el == 'G':
                return -0.5077579086
            elif el == 'H':
                return -0.5077579086
            elif el == 'I':
                return -0.5077579086
            elif el == 'K':
                return -0.5077579086
            elif el == 'L':
                return 0
            elif el == 'M':
                return -0.5077579086
            elif el == 'N':
                return -0.5077579086
            elif el == 'P':
                return -0.5077579086
            elif el == 'Q':
                return -0.5077579086
            elif el == 'R':
                return -0.5077579086
            elif el == 'S':
                return -0.5077579086
            elif el == 'T':
                return -0.5077579086
            elif el == 'V':
                return -0.5077579086
            elif el == 'W':
                return -1.015515817
            elif el == 'Y':
                return -0.5077579086
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -1.9
            else:
                while True:
                    print('Broken ' + el)
                return -100
        if place == 2:
            if el == 'A':
                return -0.6196489269
            elif el == 'C':
                return -0.6196489269
            elif el == 'D':
                return 0
            elif el == 'E':
                return -0.6601211101
            elif el == 'F':
                return -0.6196489269
            elif el == 'G':
                return -1.293791002
            elif el == 'H':
                return -0.3300605551
            elif el == 'I':
                return -0.6196489269
            elif el == 'K':
                return -0.3300605551
            elif el == 'L':
                return -0.6196489269
            elif el == 'M':
                return -0.6196489269
            elif el == 'N':
                return -0.7898698039
            elif el == 'P':
                return -0.6196489269
            elif el == 'Q':
                return -0.6196489269
            elif el == 'R':
                return -0.3300605551
            elif el == 'S':
                return 0.0
            elif el == 'T':
                return -0.3949349019
            elif el == 'V':
                return -0.6196489269
            elif el == 'W':
                return -0.6196489269
            elif el == 'Y':
                return -0.6196489269
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -1.9
            else:
                while True:
                    print('Broken ' + el)
                return -100
        if place == 3:
            if el == 'A':
                return -1.509782697
            elif el == 'C':
                return -1.509782697
            elif el == 'D':
                return -1.849974311
            elif el == 'E':
                return -1.849974311
            elif el == 'F':
                return -1.509782697
            elif el == 'G':
                return -1.509782697
            elif el == 'H':
                return -0.1937030714
            elif el == 'I':
                return -1.509782697
            elif el == 'K':
                return 0
            elif el == 'L':
                return -1.509782697
            elif el == 'M':
                return -1.509782697
            elif el == 'N':
                return -1.552458901
            elif el == 'P':
                return -1.509782697
            elif el == 'Q':
                return -1.657368255
            elif el == 'R':
                return -0.3874061428
            elif el == 'S':
                return -0.9336572298
            elif el == 'T':
                return -1.895646403
            elif el == 'V':
                return -1.509782697
            elif el == 'W':
                return -1.509782697
            elif el == 'Y':
                return -1.509782697
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -1.9
            else:
                while True:
                    print('Broken ' + el)
                return -100
        if place == 4:
            if el == 'A':
                return 0
            elif el == 'C':
                return -0.5193629289
            elif el == 'D':
                return -1.569436786
            elif el == 'E':
                return -1.058965506
            elif el == 'F':
                return -0.5193629289
            elif el == 'G':
                return -0.3374280854
            elif el == 'H':
                return -0.1487252672
            elif el == 'I':
                return -0.5193629289
            elif el == 'K':
                return -0.1487252672
            elif el == 'L':
                return -0.5193629289
            elif el == 'M':
                return -0.5193629289
            elif el == 'N':
                return -0.7396410968
            elif el == 'P':
                return -0.5193629289
            elif el == 'Q':
                return -1.220660701
            elif el == 'R':
                return -0.1487252672
            elif el == 'S':
                return -0.2131339135
            elif el == 'T':
                return -0.2886446868
            elif el == 'V':
                return -0.5193629289
            elif el == 'W':
                return -0.5193629289
            elif el == 'Y':
                return -0.5193629289
            elif el == 'Z':
                return 0.0
            elif el == '-':
                return -1.9
            else:
                while True:
                    print('Broken ' + el)
                return -100

root = tkinter.Tk()
root.withdraw()  # use to hide tkinter window

currdir = os.getcwd()

# AG names = '0002BBY,0005BBY,0008BBY,0013BBY,0014BBY,0015BBY,1475MV,1970LE,MHI,MHR,DJ263,3273,MP522,AA039,CY048,LB037,DJ258,I_2496,ITM4,KM11,KM18,RUA005,RUA010,VI1090'
# a6 names = '1016,1038,1041,1056,1080,1125,5,CY057,CY171,CY173,H386,H408,H410,K08,K84,R053,R163,R392,R497,R526,R575,R589,RUA001,RUA002,RUA003,RUA004,RUA006,RUA007,RUA009,SC1233,SC1457,SC3208,SC3410,UA0116'
# a1 names = '103,1135M,1141M,120F,120M,191084,191845,209,21,216,23,24,263,283F,320,339,366,368F,368M,374,388,391,398,405,407,426,467,503_15344_CT533,511,515,575,601F,601M,655,657,713,717,86,9004,A5340A09,BF520,CF01,CO783,Cst_013,CY153,G505mother,HCPI_13,I206mother,ITM1,KCC2,KER2008,KNH1144,MF520,ML170,ML1990,ML605,ML752,PP6_F2,PP6_F3,Q23,Q461,QA413,QB726,QB850,QF495,R463F,R774F,R880F,R880M,RW008,SF170,SK162,UG029,UG037,UG273,UG275'
# D names = '91,108F,108M,11,12,183F,183M,190049,191647,191727,191859,191882,2769_F,2769_M,2810_F,2810_M,295F,295M,326F,326M,32F,32M,338F,338M,372F,372M,385,394_F,394_M,527,54,552,57128,602F,602M,605F,605M,888_F,888_M,890_F,890_M,927_F,927_M,A03349,A03836_9F,A07412,A08483,B03887_9M,D053826,DM1,DM2,DM3,DM4,DM5,DM6,DM7,DM8,DM9,F053828,HF_F,HF_M,J32228,KBH,P4b,QA013,QA465,QD435,RJ100,SE365,SR20_F,SR20_M,SR5_F,SR5_M,UG001,UG024,UG046,UG065,UG270,UG274,ZM387'

# B Eur names = '131,133,153,159,243,309,19298,19554,19642,19663,19956,107,111,1HD10K,1HD11K,1HD4K,1HD5K,1HD6K,1HD8K,1HD9K,A,ACH18969,ACH19542,ACH19659,ACH19768,ACH19999,B,C,CW002,CW010,CW012,CW048,D,E,E21,LAI,LW,MDM,MM4,MM8,N2,NA20,NAB01,P1,patient,SHCS566,SHCS604,SHCS617,SHCS901,UKBH1,ZPHI11,ZPHI127,ZPHI143,ZPHI149,ZPHI158,ZPHI24,ZPHI39,ZPHI67,ZPHI72,ZPHI92'
#AE names = '3002,3019,3022,3025,3043,3046,3063,3067,3078,3084,3090,3104,3111,3112,3118,3131,3135,3151,3152,3153,3162,3184,3186,3189,3193,3203,3210,3212,3218,3219,254002,254004,254016,254019,254023,254024,254029,254032,254036,320041,AA001,AA002,AA003,AA004,AA005,AA006,AA007,AA008,AA009,AA012,AA013,AA014,AA016,AA017,AA018,AA019,AA021,AA022,AA023,AA024,AA026,AA027,AA028,AA029,AA030,AA031,AA032,AA033,AA035,AA036,AA037,AA041,AA042,AA044,AA045,AA048,AA049,AA050,AA051,AA052,AA054,AA057,AA058,AA059,AA061,AA062,AA063,AA064,AA065,AA066,AA067,AA068,AA069,AA070,AA072,AA073,AA074,AA075,AA076,AA077,AA078,AA079,AA080,AA081,AA082,AA083,AA086,AA088,AA089,AA090,AA091,AA092,AA094,AA097,AA098,AA099,AA100,AA101,AA102,AA103,AA104,AA105,AA107,AA108,AA109,AA110,AA111,AA112,AA113,AA116,AA117,AA118,AA119,AA121,AA122,AA123,AA125,AA126,AA127,AA129,AA130,AA131,CA10,CM235,CM240,CM244,CR02,CR03,CR08,CR10,CR11,CR14,CR15,CR19,CR25,CR28,CR36,CR38,ID17,NH1,OX005000,OX008000,OX009000,OX010000,OX012000,OX015000,OX017000,OX018000,OX021000,OX023000,OX025000,OX028000,OX030000,OX031000,T500104_256254,T500105_293735,T500107_357545,T500107_535902,T500204_502281,T500206_614109,T500207_502102,T500207_503006,T500207_509989,T500207_518020,T500208_504258,TH001,TH011'
#C AFrica names='89,114,334,393,478,595,626,665,682,985,1086,1172,1176,1196,1335,1373,1394,2010,2052,2060,2103,3002,3003,3004,3006,3009,3017,3032,3036,3037,3040,4001,4002,4004,4007,4008,4013,4014,4015,4016,4017,4026,4027,4028,4029,4030,4031,4032,4034,4036,4037,4039,4041,4045,4046,4048,4049,4050,4055,4056,4058,4059,4061,702010141,703010010,703010054,703010131,703010159,703010167,703010193,703010200,703010217,703010228,703010269,714903009,714903305,714903404,714903902,714904707,714910404,715900301,715900702,C007,C009,C011,C012,C018,C019,C030,C034,C047,C059,C061,C070,C083,C109,QB008,QD022,SM145'
#names = '1,18,21,22,23,25,26,28,135,148,154,1001,1006,1012,1018,1053,1058,1059,1674,1675,12007,12008,1HB1,1HB2,1HB3,1HC1,1HC2,1HC3,1HD1,2A1,2C5,2e1,2e2,2e4,2e5'
#B Cari names = 'SC02,SC05,SC11,SC13,SC20,SC22,SC25,SC31,SC33,SC42,SC45,SC46,SC50,TT103,TT106,TT112P,TT113P,TT114P,TT27P,TT28P,TT29P,TT31P,TT34P,TT35P'
#USCC names = '1,18,21,22,23,25,26,28,135,148,154,1001,1006,1012,1018,1053,1058,1059,1674,1675,12007,12008,1HB1,1HB2,1HB3,1HC1,1HC2,1HC3,1HD1,2A1,2C5,2e1,2e2,2e4,2e5,SC02,SC05,SC11,SC13,SC20,SC22,SC25,SC31,SC33,SC42,SC45,SC46,SC50,TT103,TT106,TT112P,TT113P,TT114P,TT27P,TT28P,TT29P,TT31P,TT34P,TT35P'
#HCV _with hyphen Acute names = '6213-1,6222-08,10002-04,10002-07,10012-06,10012-08,10012-1,10012-13,10017-09,10017-1,10017-12,10017-14,10017-16,10020-08,10020-13,10021-1,10021-14,10024-06,10024-07,10024-08,10025-08,10025-09,10025-11,10029-08,10029-09,10029-11,10029-15,10051-11,10051-14,10062-03,10062-04,10062-05,10062-08'
#HCV Chronic names = 'ARJA6267,BLMI6862,BRRO6924,GOTO90002,JOTO6422,KNPH3730,LAST90001,ROMI6847,RUVI5913,SLRO5563,WEPA5774,WHRO3882,WIMI4025,WIMI90003'
#HCV Acute names = '6213,6222,10002,10012,10017,10020,10021,10024,10025,10029,10051,10062'
#HCV clade3 names = '9055,10003,ILBSRAS,K710,patient,Pt4,Pt5,R709,RASILBS,S48,S49,UKN3A1,UKN3A2,UKN3A4'
#Longitudinal AE names= '3022,3025,3046,3090,3112'
#Longitudinal B names = '1,133,135,159,169,19298,19554,19642,19663,19956,2008,2016,2019,25,26,28,3002,309,4048,45,700010040,700010058,700010077,700010470,7115,9,9006,9018,9021,9024,9040,9213,ACH142,ACH18969,ACH19542,ACH19659,ACH19768,ACH19999,ES2,ES8,ES9,H434,H671,HIA,HP11,HP14,HP16,HP17,HP18,HP19,HP2,HP3,HP4,HP5,HP7,HP8,HP9,IIA,KYR,LIA,LW,MM42,MM43,mother,P1024,P1031,P1046,P1189,patient,PIC1362,PIC55751,PIC83747,PIC90770,POS1,POS2,SC02,SC05,SC11,SC13,SC24,SC25,SC51,VC10014,WC1,WC3,WEAU0575,Z258'
#B_NA2 names = '306159,306344,306376,306512,306517,4012,4013,4013171,4013211,4013240,4013242,4013291,4013296,4013321,4013327,4013383,4013396,4013419,4013440,4013446,4013448,4030,4033,4051,4059,5002,5003,502-0879,546,546,57,61,61792,61792,62130,62130,62357,62357,6240,6244,6247,6248,62615,62615,62995,62995,63054,63054,63068,63068,63215,63215,63358,63358,63396,63396,6535-P30,6568,6568,700010019,700010025'

#AllLong USCC
#names = '1,9,25,26,28,45,135,169,700010040,700010058,700010077,700010470,ES2,ES8,ES9,mother,P1024,P1031,P1046,P1189,PIC1362,PIC55751,PIC83747,PIC90770,POS1,POS2,SC02,SC05,SC11,SC13,SC24,SC25,SC51,VC10014,WC1,WC3,WEAU0575,Z258,1051,1058,1059,148,22,23,6247,700010148,9001,B106,B115,B155,B199,BORI0637,SUMA0874'
#Long Europe
#names = '133,159,19298,19554,19642,19663,19956,309,9213,ACH142,ACH18969,ACH19542,ACH19659,ACH19768,ACH19999,H434,H671,LW,MM42,MM43,patient'


#NEW DATA

#AE_Acute
#names ='LAI,3002,3019,3022,3025,3043,3046,3063,3067,3078,3084,3090,3104,3111,3112,3118,3131,3135,3151,3152,3153,3162,3184,3186,3189,3193,320041,3203,3210,3212,3218,3219,AA001,AA002,AA003,AA004,AA005,AA006,AA007,AA008,AA009,AA012,AA013,AA014,AA015,AA016,AA017,AA018,AA019,AA021,AA022,AA023,AA024,AA026,AA027,AA028,AA029,AA030,AA031,AA032,AA033,AA034,AA035,AA036,AA037,AA038,AA041,AA042,AA044,AA045,AA048,AA049,AA050,AA051,AA052,AA054,AA055,AA057,AA058,AA059,AA060,AA061,AA062,AA063,AA064,AA065,AA066,AA067,AA068,AA069,AA070,AA072,AA073,AA074,AA075,AA076,AA077,AA078,AA079,AA080,AA081,AA082,AA083,AA085,AA086,AA088,AA089,AA090,AA091,AA092,AA094,AA097,AA098,AA099,AA100,AA101,AA102,AA103,AA104,AA105,AA107,AA108,AA109,AA110,AA111,AA112,AA113,AA116,AA117,AA118,AA119,AA121,AA122,AA123,AA125,AA126,AA127,AA129,AA130,AA131,OX005000,OX008000,OX009000,OX010000,OX012000,OX015000,OX017000,OX018000,OX021000,OX023000,OX025000,OX028000,OX030000,OX031000'

#AE Chronic
#names ='CMU10,CR10,NH1,T500104_256254,T500105_293735,T500107_357545,T500107_535902,T500204_502281,T500206_614109,T500207_502102,T500207_503006,T500207_509989,T500207_518020,T500208_504258'

#B_NA_Acute
#names = 'LAI,1001,1006,1012,1018,1053,12007,12008,4013171,4013211,4013240,4013242,4013291,4013296,4013321,4013327,4013383,4013396,4013419,4013440,4013446,4013448,61792,62130,62357,6240,6244,6248,62615,62995,63054,63068,63215,63358,63396,700010081,700010106,700010224,700010238,700010246,700010607,700010610,700010649,700010685,700010867,700010914,700010937,701010016,701010027,701010043,701010055,701010068,701010092,701010302,701010344,701010378,701010449,9010,9014,AD17,AD18,AD75,AD77,AD83,CAAN5342,INME,MEMI4948,PRB931,PRB956,PRB957,PRB958,PRB959,REJO4541,RHPA4259,THRO4156,TRJO4551,US006,WITO4160,Z02,Z03,Z05,Z10,Z13,Z16,Z18,Z20,Z23,Z27,Z29,Z30,Z31,Z32,Z34,Z35,Z62,Z64,Z71,Z74,Z75,Z78,Z85,Z86,Z91,Z92,Z93,Z94,Z95,700010654,AS1,N8261,T8250,306159,306344,306376,306512,306517,9002,9003,9016,9017,9022,9023,9025,9026,9027,9028,9029,9030,9031,9032,9033,9036,9044,9045,9048,9055,9058,9062,9063,9071,9075,9076,9077,9079,9082,9083,9096,9097'


#C_Eastern_Central_Acute_1
#names = 'LAI,114,334,393,478,595,626,665,682,89,QB008,QD022,703010193,703010200,703010217,703010228,703010010,703010054,2899_Partner,2173_Partner,3017,3037'


#C_Eastern_Central_Acute_2
#names = 'LAI,14B,17B,18B,1B,234,304,390,3B,410,556,569,5B,6B,7B,9B,Z1792M,ZM1072M,ZM1464M,ZM1781M,ZM178F,ZM180M,ZM184F,ZM197M,ZM206F,ZM214M,ZM215F,ZM229M,ZM233M,ZM235F,ZM237M,ZM247M,ZM267F,ZM284M,ZM289M,ZM297M,ZM503F,ZM246M'



#C_Southern_Africa_Acute
#names = 'LAI,1245045,19157834,20258279,20927783,21197826,21283649,2833264,2935054,503_02660_CT184,503_06150_CT885,503_06877_KO870,503_08252_SO405,503_09003_CA392,503_13503_SO722,503_13580_CA146,503_15383_Me178,704010017,704010056,704010069,704010083,705010015,705010026,705010078,705010110,B001811,B005018,B005582,CAP129,CAP136,CAP174,CAP188,CAP200,CAP217-dec2005,CAP221,CAP222,CAP224,CAP237,CAP63,CAP69,CAP84,Du36,704809221,704810053,706010018'


#B_NA_Chronic
#names = '4013226,1018,1674,1675,1775,18,1825,1853,4012,4013,4030,4033,4051,4059,5002,5003,502-0879,57,61,700010025,700010135,700010150,700010174,700010252,700010260,700010287,700010464,700010470,700010728,700010734,700010742,700010791,700010920,701010114,701010504,701010533,701010541,7036,85,A8110,AVDA2874,BEKA5842,BRJO4843,C109,C62,C93,C94,C96,C98,COCE6096,CRPE4571,DIGA3757,EABE4469,ES10,ES3,ES4,ES5,ES6,ES7,F6817,FERI5029,FOJO4081,GETO5098,HALA6323,HEMA4284,J7180,JACH1853,JOSA5789,JOTO5278,K8072,LAHA4867,MCBR4209,MCRO3633,MCST4474,MEJA5586,OLLA4645,RHGA1581,RHMI4089,RIER,ROCH4447,ROST4216,SADO6038,SAMI4303,SHKE4761,SMRE4166,SPFE4120,STCO5453,T4590,T8107,TALA4022,UNC2009,UNC3405,UNC4295,UNC4484,UNC4911,UNC5057,UNC5283,UNC5417,UNC5479,UNC5539,UNC5548,UNC5734,UNC5769,UNC5791,UNC5799,UNC6064,UNC7092,UW1218,UW1249,UW1330,UW1352,UW1423,UW1444,UW1446,UW1451,UW1470,UW1508,UW1586,UW1588,UW1599,UW1624,UW1631,UW1632,UW1711,UW1791,UW1794,WARO5662,WICU4248,YOMI4024,9076'

#B_Europe_Chronic
#names = '20044616,CW002,CW010,CW012,CW048,MM4,MM8,NL43,Pat-1,PVO,SHCS396,SHCS566,SHCS604,SHCS617,SHCS901,ZPHI18,ZPHI67'

#C_Eastern_Central_Chronic_1
#names = '1168_Index,1168_Partner,2173_Index,2899_Index,3003,3004,3011,3012,3022,3025,3026,3027,3029,3031,3032,3034,3036,3040,3042,3045,3048,3050,702010141,703010167,703010269,714903009,714903305,714903404,714904707,714905807,714910404,715900301,715900702,715901209,C111,C113,C120'

#C_Eastern_Central_Chronic_2
#names = '053M,055F,10M,11M,12M,135M,14M,16M,17M,19M,31M,32M,33M,3M,4M,5M,6M,7M,9M,ZM375,ZM376,ZM377,ZM378,ZM379,ZM381,ZM382,ZM383,ZM384,ZM388,ZM389,ZM393,ZM394,ZM395,ZM399,ZM400,ZM401,ZM402,ZM403,ZM405,ZM406,ZM408,ZM410,ZM411,ZM412,ZM413,ZM414,ZM415,ZM416,ZM418,ZM419,ZM420'


#C_Southern_Africa_Chronic
#names = '704010028,704010034,704010075,704010207,704010273,704010330,704010499,705010110,CAP88,CF04,CF05,CF08,CF09,CF10,CF13,CF17,CF21,Du10,DU114,DU179,Du246'

#IC_Chronic
#names = 'IC_10048,IC_10196,IC_10204,IC_10213,IC_10259,IC_10425,IC_30048,IC_20476,IC_10419,IC_10473,IC_11100,IC_11200,IC_20439,IC_30281,IC_27,IC_38,IC_61,IC_101,IC_102,IC_122,IC_123,IC_18,IC_193,IC_2,IC_23,IC_231,IC_239,IC_241,IC_246,IC_247,IC_253,IC_54,IC_254,IC_53,IC_35,IC_387,IC_389,IC_404,IC_434,IC_450,IC_451,IC_454,IC_461,IC_462,IC_465,IC_479,IC_277,IC_357,IC_525,IC_616,IC_619,IC_645,IC_689,IC_714,IC_715,IC_295,IC_537,IC_736,IC_741,IC_742,IC_691,IC_818,IC_206,IC_665,IC_684,IC_711,IC_786,IC_797,IC_912,IC_935,IC_940,IC_1036,IC_660,IC_999,IC_25,IC_449,IC_480,IC_534,IC_936,IC_1013,IC_629,IC_798,IC_1115,IC_1211,IC_446,IC_857,IC_1125,IC_1170,IC_1240,IC_1260,IC_493,IC_606,IC_907,IC_1329,IC_1331,IC_1336,IC_1378,IC_1386,IC_472,IC_520,IC_DAH,IC_JBS,IC_JCP,IC_JJM,IC_REH'

#UW_Chronic
#names = 'UW_1002,UW_1295,UW_1375,UW_1456,UW_1555,UW_1641'



#Longitudinal names
#B_NA_LS_4
#names = 'LAI,26,700010040,700010058,700010077,700010470,ES8,ES9,PIC1362,PIC55751,PIC83747,PIC90770,SC05,SC13,SC24,SC51,WC3,WEAU0575'

#B_NA_LS 2 3
#names = 'LAI,1051,1058,1059,6247,9001,BORI0637,SUMA0874,VC10014'

#IC_LS
#names = 'IC_473,IC_656,IC_790,IC_475,IC_510,IC_669,IC_663'

#UW_LS
#names = 'UW_1140,UW_1261,UW_1266,UW_1313,UW_1368,UW_1386,UW_1393,UW_1406,UW_1842'

#B_Eur_LS
#names = 'LAI,133,159,309,9213,LW,patient'

#C_SA_LS
#names = '704010042,705010162,705010185,705010198,706010164,CAP177,CAP210,CAP239,CAP256,CAP45,Du123,DU151'

#C_ECA_LS
#names = '221M,703010131,703010159,703010256,703010505,707010457,ZM1024F,ZM153M,ZM247F'



#Treated

#B_Eur_RX:
#names = '1HD10K,1HD11K,1HD4K,1HD5K,1HD6K,1HD8K,1HD9K'

#B_NA1_RX:
#names = '21,2A1,2C4,2C5,2E1,2E2,2E4,2E5'

#B_NA4_RX:
#names = 'LAI,A5340A01,A5340A02,A5340A03,A5340A04,A5340A05,A5340A06,A5340A07,A5340A08,A5340A09,A5340A10'

#Temp for Rentian
names = '1HB1,1HB2,1HB3,1HC1,1HC2,1HC3,1HD10K,1HD11K,1HD1,1HD2,1HD4K,1HD5K,1HD6K,1HD8K,1HD9K'


name = names.split(',')
for nam in name:
    FiveD = False
    weight_scale = 'G' #F for 2F5 or G for 2G12
    try:
        #distances = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Volatility_Data\\V3 Loop Volatility\\For Volatility- Myside\\Distances\\HCV\\Chronic\\" + nam + '.fas')
        #distances = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Analysis\\New Project\\Genetic Distances\\B_Europe\\" + nam + '.fas');
        distances = open(
            "C:\\Users\\kandula.HEALTHCARE\\Desktop\\Temp\\TempForRentian\\Long_2_B_NA1\\" + nam + '.fas');

    except:
        print(nam)
        continue
    #table = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Analysis\\Hepatitis C\\Volatility\\Data\\Tab Delimited\\NewPositions\\Chronic.txt")
    #D file
    #table = open("U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Analysis\\New Project\\Volatility Results\\B_LS_Europe_ForVolatility.txt")
    table = open("C:\\Users\\kandula.HEALTHCARE\\Desktop\\Temp\\TempForRentian\\tempFile.txt")
        #"U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Analysis\\New Project\\AllClades_JustSequences\\Glyco Sites_All Clades\\TSV\\UW_Chronic.txt")
        #"U:\\ResearchData\\rdss_hhaim\\LAB PROJECTS\\Raghav\\Analysis\\Env Volatility Forecasting Project\\EVF MS Data and Analyses\\Working Directory\\TSV files_PNGS\\C_Southern_Africa_Chronic.txt

    # HyPhy Philips Style Matrix has number of Sequences listed in the first position of output#
    # This code brings the matrix into python#
    count = 0
    tot = {''}
    tot.remove('')
    dm = {}

    for row in distances:
        if count == 0:
            count += 1
            continue
        t = row.split()
        if nam[:2] == '2e':
            t[0] = t[0][2:]
            t[1] = t[1][2:]
            wait1 = False
            fixer = ''
            for i in t[0]:
                if i == '-':
                    wait1 = True
                elif wait1:
                    if i == '.':
                        wait1=False
                        fixer += '.'
                else:
                    fixer += i
            z0=fixer
            wait1 = False
            fixer = ''
            for i in t[1]:
                if i == '-':
                    wait1 = True
                elif wait1:
                    if i == '.':
                        wait1=False
                        fixer+='.'
                else:
                    fixer += i
            z1 = fixer
        else:
            fixer = ''
            for i in t[0]:
                if i == '-':
                    pass
                else:
                    fixer += i
            z0=fixer

            fixer = ''
            for i in t[1]:
                if i == '-':
                    pass
                else:
                    fixer += i
            z1 = fixer


        if z0 not in tot:
            tot.add(z0)
        if z1 not in tot:
            tot.add(z1)
        dm[(z0, z1)] = t[2]
        dm[(z1, z0)] = t[2]
        count += 1
    numRows = (len(tot))
    if numRows < 2:
        continue
    #################################################################################


    weights = np.array([1, 1, 1, 1, 1])
    #weights = np.array([1.237176432,0.6845558143,0.2308496402,0.8544987065,0.7873302248]) #2F5
    #weights = np.array([1.981905578,2.017661388,0.712722561,2.232499063,0.5831526609]) #2G12



    # Brings table into Python #
    lines = [[]] * (numRows + 1)
    count = 0

    years = []
    phases = []
    for line in table:
        if len(line) == 0:
            break
        t = line.split()
        #print(line);
        if t[0] == 'Country':
            lines[count] = t
            count += 1
            length = len(t) - 4 #number of items in the name_space
            continue

        name = t[0] + '.' + t[1] + '.' + t[2].split("-")[0] + '.' + t[3]
        if name not in tot:
            #print(name)
            continue
        if t[1] not in years:
            years.append(t[1])
        try:
            if t[2].split("-")[1] not in phases:
                phases.append(t[2].split("-")[1])
        except:
            if "" not in phases:
                phases.append("")
        lines[count] = t
        count += 1
    Positions = [0] * (length)

    c = 0
    for pos in lines[0]:
        if len(pos) != 3:
            continue
        Positions[c] = pos
        c += 1

    count = 0

    num = 0
    vols = [0] * len(Positions)
    cont = False
    for year in years:
        if year == 'None':
            continue
        if cont:
            break
        for phase in phases:
            posCount = 0
            for pos in Positions:
                EnvIDs = [0] * numRows
                PatNum = [0] * numRows
                positionVals = [0] * (numRows)
                count = 0
                yr_counter = 0
                for line in lines:
                    skipper = False
                    already_skip = True
                    if len(line) == 0:
                        break
                    if line[1] != year:
                        yr_counter += 1
                        skipper = True
                        already_skip = False
                    try:
                        if line[2].split("-")[1] != phase:
                            if already_skip:
                                yr_counter += 1
                                skipper = True
                    except:
                        ZZ = 0

                    if count < 1:
                        count += 1
                        continue
                    EnvIDs[count - 1] = line[0] + '.' + line[1] + '.' + line[2].split("-")[0] + '.' + line[3]
                    if(len(line) > (posCount+4)):
                        #print(line)
                        #print(nam)
                        #print(posCount)
                        #print(len(line))
                        #Added this condition - Raghav
                        positionVals[count - 1] = getVal(line[posCount + 4]) / 100 #hydropathy
                        #positionVals[count - 1] = getVal(line[posCount + 4]) / 100  # hydropathy
                        #positionVals[count - 1] = vals(weight_scale,posCount, line[posCount+4]) #weight-scaling
                    if FiveD:
                        positionVals[count - 1] = np.array([getVal(line[posCount + 4]) / 100, getVal(line[posCount + 5]) / 100, getVal(line[posCount + 6]) / 100, getVal(line[posCount + 7]) / 100, getVal(line[posCount + 8]) / 100])
                    if skipper:
                        if FiveD:
                            positionVals[count - 1] = np.array([-1000])
                        else:
                            positionVals[count - 1] = -1000
                    count += 1

                deltaVals = [[0.0 for x in range(numRows)] for y in range(numRows)]
                volVals = [[0.0 for x in range(numRows)] for y in range(numRows)]
                totalVol = 0
                numCalcs = 1
                for i in range(0, numRows):
                    for j in range(0, numRows):
                        if i == j:
                            break
                        if FiveD:
                            try:
                                if positionVals[i][0] == -1000 or positionVals[j][0] == -1000:
                                    continue
                            except:
                                pass
                        else:
                            if positionVals[i] == -1000 or positionVals[j] == -1000:
                                continue
                        diffs = (positionVals[i] - positionVals[j])
                        if FiveD:
                            weighted = LA.norm(diffs)
                            deltaVals[i][j] = weighted
                        else:
                            deltaVals[i][j] = math.pow(diffs, 2)
                        try:
                            t = float(dm[(EnvIDs[i], EnvIDs[j])])
                            r = round(float(dm[(EnvIDs[i], EnvIDs[j])]), 8)
                            if r > 0:
                                volVals[i][j] = round(deltaVals[i][j] / float(dm[(EnvIDs[i], EnvIDs[j])]), 8)
                                t = round(deltaVals[i][j] / float(dm[(EnvIDs[i], EnvIDs[j])]), 8)
                                totalVol += volVals[i][j]
                                numCalcs += 1
                        except:
                            pass

                vols[posCount] = str(totalVol / numCalcs)
                posCount += 1

                num = count - yr_counter

                if FiveD:
                    break
            # if cont:
            #     continue
            vol_string = " "
            for vol in vols:
                vol_string += str(vol) + ' '

            try:
                if num > 1:
                    if FiveD:
                        print(str(num) + ' ' + line[0] + ' ' + year + ' ' + line[2].split("-")[0] + ' ' + phase + ' ' + line[3] + '  ' + vols[0])
                    else:
                        print(str(num) + ' ' + line[0] + ' ' + year + ' ' + line[2].split("-")[0] + ' ' + phase + ' ' + line[3] + ' ' + vol_string)
            except:
                if FiveD:
                    print(str(num) + " " + '-' + ' ' + year + ' ' + nam + ' ' + phase + ' ' + '-' + '  ' + vols[0])
                else:
                    print(str(num) +" "+ '-' + ' ' + year +' '+ nam +' '+ phase+' '+'-'+' '+vol_string)