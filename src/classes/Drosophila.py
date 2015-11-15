'''
Created on 16/05/2013

@author: Flavio Lichtenstein
@local: Unifesp - Bioinformatica

'''
from Bio import Entrez

class Drosophila:
    def __init__(self):
        self.mnem =  ['amer', 'amer2', 'amertex', 'ana', 'ariz', 'kik', 'mau', 'mel', 'mir', 'moja', 'pseud', 'pbog', 'pseud2', 'sec', 'sim', 'subob', 'teis', 'viri', 'wil', 'yak']

        self.speciesList = ['paulistorum', 'willistoni',\
    'pseudoobscura bogotana', 'pseudoobscura pseudoobscura', 'pseudoobscura',\
    'sinobscura','subobscura', 'obscuripes', 'obscura',\
    'pseudoananassae','prolongata','curveadeagus', 'kohkoa', 'annulipes', \
    'ohnishii','mayri','liui','ambigua', 'oritisa', 'curviceps', 'bunnanda', \
    'bocki', 'baimaii','serrata', 'malerkotliana', 'ornatifrons', 'mesophragmatica', 'brncici', \
    'kuntzei', 'madikerii','takahashii', 'ogumai','ficusphila', 'gasici', 'viracochi', \
    'trapezifrons', 'prostipennis', 'ornatipennis', 'curveadeagus','bipectinata intronless', \
    'pulchrella','watanabei', 'trilutea', 'tani', 'subauraria',  \
    'suzukii','seguyi','ficusphila','yakuba','subelegans', 'gaucha', \
    'madeirensis','affinidisjuncta','erecta','sechellia','mojavensis','simulans',\
    'persimilis','virilis','grimshawi','ananassae','vallismaia','takahashii',\
    'palustris', 'latifshahi', 'hydeoides', 'daruma', 'angularis', 'adamsi', 'ustulata',\
    'subpalustris', 'mimica', 'kanekoi', 'aracea','fima','eugracilis','bipectinata',\
    'aracataca', 'watanabei', 'suzukii', 'seguyi','pseudoananassae','prolongata','ohnishii',\
    'mayri', 'liui', 'curveadeagus', 'bocki', 'baimaii', 'malerkotliana', 'serrata',\
    'madikerii', 'ogumai', 'constricta', 'lutescens', 'tenebrosa', 'limingi', \
    'dianensis', 'luguensis', 'hubeiensis', 'hypercephala', 'mauritiana', \
    'tetravittata', 'littoralis', 'mercatorum', 'testacea',\
    'busckii', 'hydei', 'limensis', 'mercatorum', 'nigrodumosa', 'bifurca', \
    'parvula', 'fuyamai', 'birchii', 'limingi', 'polymorpha', 'phaeopleura', 'pallidipennis',\
    'ochrogaster', 'nigromaculata', 'neomorpha', 'nebulosa', 'mimetica',\
    'melanogaster', 'nigromelanica', 'melanocephala', 'melanica',\
    'caribiana', 'ararama', 'lummei', 'kepulauana', 'sturtevanti', 'siamana', 'rubida', 'wheeleri',\
    'sulfurigaster', 'polychaeta', 'phalerata', 'pavani', 'pallidifrons', 'novamexicana',\
    'neocordata', 'lachaisei', 'nasuta', 'merina', 'mediopictoides', 'flavohirta', 'bromeliae',\
    'albirostris', 'hypocausta', 'gibberosa', 'orena', 'varians','algonquin',\
    'repletoides', 'bicornuta', 'quadraria', 'pallidosa', 'punjabiensis', 'immigrans',\
    'levii', 'greeni', 'camargoi', 'albomicans', 'biarmipes', 'formosana', 'cardini', 'lucipennis',\
    'punjabiensis', 'nikananu', 'leontia', 'chauvacae', 'nagarholensis', 'monieri', 'malagassya',\
    'funebris', 'arawakana', 'repleta', 'rufa', 'biauraria', 'kitumensis', 'helvetica','white tip scutellum',\
    'subsilvestris', 'silvestris', 'planitibia', 'microlabis', 'tristis', 'bifasciata', 'affinis', 'athabasca','soonae',\
    'imaii', 'bakoue', 'atripex', 'auraria', 'jambulina', 'ercepeae', 'elegans','dunni', 'narragansett',\
    'asahinai', 'davidi', 'cauverii', 'diplacantha', 'tropicalis', 'tsacasi', 'vulcana', 'miranda','neokadai',\
    'metzii','adunca','tanythrix','petalopeza','lineosetae','dasycnemia','cyrtoloma','borealis','angor',\
    'tolteca',  'azteca', 'narragansett', 'tsukubaensis', 'dossoui', 'lini', 'kikkawai','sp. T',\
    'bocqueti', 'sternopleuralis', 'santomea', 'ruberrima', 'tsigana', 'transversa', 'talamancana','disjuncta','velox',\
    'quadrisetata', 'fluvialis', 'sucinea', 'insularis', 'equinoxialis', 'pavlovskiana', 'fumipennis', 'capricorni',\
    'americana texana', 'americana americana', 'americana', 'texana', 'flavomontana', 'montana',\
    'differens', 'hawaiiensis', 'lacicola', 'nigra', 'adiastola', 'mettleri', 'picticornis',   \
    'heteroneura', 'karakasa', 'beppui', 'perlucida', 'asper', 'cf. cheda', 'nannoptera', 'ezoana', 'anceps', \
    'lacertosa', 'yunnanensis', 'bipatita', 'multidentata', 'pseudosordidula', 'sordidula', 'moriwakii', \
    'unimaculata', 'robusta', 'okadai', 'euronotus', 'emarginata', 'milleri',  'pinicola', \
    'subsaltans', 'austrosaltans', 'lusaltans','prosaltans','saltans', 'unispina', 'curvispina', \
    'orosa', 'kanapiae', 'paralutea', 'tripunctata', 'ingens', 'setosifrons', 'substenoptera', 'neoperkinsi', \
    'nigribasis', 'nigrilineata', 'primaeva', 'gunungcola', 'quadrilineata', 'maculinotata', 'guttifera', \
    'endobranchia', 'antillea', 'hanaulae', 'hemipeza', 'oahuensis', 'neopicta', \
    'papuensis',\
    'limbata','teissieri','guanche','iri', 'iki','guaru', 'hei', 'burlai', 'gani', 'eskoi', 'barbarae', 'bai']

        ''' list of possible Gene names '''
        self.adhList = ['ADH-1', 'ADH1', 'ADH-2', 'ADH2', 'ADH-F', 'ADHF', 'ADH-S', 'ADHS', 'ADH-R', 'ADHR', 'ADH']

        self.amyList = ['AMY-P', 'AMY-D', 'AMYREL', 'AMY']


        self.listSpeciesLittle = ['americana', 'americana americana', 'ananassae', 'arizonae', 'borborema', 'gouveai', 'innubila', 'kikkawai', 'mauritiana', 'melanogaster',  'miranda', 'mojavensis', 'pachea', 'pseudoobscura',  'santomea', 'sechellia', 'simulans', 'subobscura', 'teissieri', 'virilis', 'yakuba', 'willistoni' ]


        self.listGenes = ['128up', '140up', '17.6env', '17.6gag', '18w', '297pol', 'A3-3', 'a6', 'Aats-glupro', 'aay', 'abd-A', 'abo',  \
'abs', 'Ace', 'acj6', 'Acp1', 'Acp134', 'Acp157a', 'Acp158', 'Acp16a', 'Acp16b', 'Acp16c', 'Acp19', 'Acp2', \
'Acp21a', 'Acp223', 'Acp225', 'Acp23D4a', 'Acp23D4b', 'Acp23D4c', 'Acp24', 'Acp24A4', 'Acp25', 'Acp26Aa', \
'Acp26Ab', 'Acp27a', 'Acp27b', 'Acp29AB', 'Acp29ab', 'Acp3', 'Acp32CD', 'Acp33A', 'Acp36DE', 'Acp42',  \
'Acp48', 'Acp53C14a','Acp53C14b', 'Acp53C14c', 'Acp53Ea', 'Acp54A1', 'Acp5a', 'Acp5b', 'Acp62F', 'Acp63F', \
'Acp7', 'Acp70A', 'Acp70Aa', 'Acp70Ab', 'Acp76A','Acp8', 'Acp95EF', 'Acp98AB', 'Acph', 'Acph-1', 'actinin', \
'Actn', 'Acyp', 'ade3', 'ade5', 'Adh', 'Adh-1', 'Adh-2', 'Adh1', 'Adh2', 'AdhPsi', 'Adhr', 'AdhT', 'Ag5r',  \
'AGO1', 'agt', 'al', 'Alas', 'Ald','ald','Aldh','Alk','alpha-Est1','alpha-Est10','alpha-Est2','alpha-Est3', \
'alpha-Est4','alpha-Est5','alpha-Est6','alpha-Man-I','alpha-Man-II','alpha-Spec','alphaTry','AlstR','aly', \
'Ama','amd','Amy','Amy-d','Amy-p','amy-p','Amy1','Amy2','Amy4N','Amyc1','Amyc6','Amyrel','angel','AnnX', \
'anon-25D26Ab','anon-EST:fe1G5','Anp','Antp','Anxb11','Aos1','ap','AP-50','ApepP','Appl','Arf102F','Arf84F', \
'Arp66B','Arp87C','Art1','asf1','ash1','asp','Ast','Asx','Atg9','ato','AttA','AttC','Atu','aub','aur','Axn', \
'Axs','B-H1','B-H2','bab2','bai','bam','Bap','bap','Bap55','baz','bcd','BCL7-like','bcn92','Bem46','ben',\
'beta-Spec','betaggt-I','betaTry','betaTub60D',"beta'Cop",'BG4','bgcn','bi','bib','bip1','bip2','Bj1','bnb', \
'BobA','bor','boss','br','brat','brca2','Brd','Brd8','brn','Bruce','bsk','bt','btl','btn','Bub3','burs','bw', \
'Bx','Bx42','Bzd','cac','cact','Cad96Ca','Caf1','CanA-14F','CanA1','CanB2','caps','cato','Catsup','cav','CBP', \
'Cbp53E','CCKLR-17D1','Ccp84Aa','Ccp84Ab','Ccp84Ac','Ccp84Ad','Ccp84Ae','Ccp84Af','Ccp84Ag','CDC2','cdc2', \
'Cdc37','CDC42','CDC45L','cdc6','Cdic','Cec-Psi1','Cec1','Cec2','CecA2','CecB','CecC','Cen190','CG10035', \
'CG10198','CG10200','CG10252','CG10260','CG10307','CG10321','CG10470','CG10521-like','CG10561','CG10623', \
'CG10680','CG10750','CG10853','CG10912','CG10920','CG10990','CG10990-like','CG10996','CG11023','CG11037', \
'CG11092','CG11093','CG11105','CG11105-like','CG11126','CG11152','CG11190','CG11354-like','CG11378','CG11379', \
'CG11380','CG11382','CG11384','CG11387-like','CG11398','CG11403','CG11426','CG11473','CG11475','CG11638', \
'CG11664','CG11697','CG11715-like','CG11727','CG11779','CG11802','CG11802-like','CG11971','CG12123','CG12139', \
'CG1217-like','CG12199','CG12244-like','CG12314','CG1239','CG12395','CG12470','CG12496','CG12498','CG12608', \
'CG12681','CG12720','CG12737-like','CG12772','CG12780','CG12909','CG13004','CG1307','CG1314','CG1315','CG13189', \
'CG13277','CG13350','CG13375','CG13377','CG13422','CG13527','CG13617','CG13690','CG13732','CG13758','CG13759', \
'CG13760','CG13762','CG13845','CG13895','CG13934','CG1397','CG14045','CG14045-like','CG14050','CG1440', \
'CG14408','CG14414','CG14416','CG14417','CG14418','CG14434','CG14435','CG14435-like','CG14438','CG14624', \
'CG14625','CG14626','CG14627','CG14628','CG14629','CG14630','CG14717','CG14785','CG14786','CG14787','CG14797', \
'CG14806','CG14926','CG14957','CG1503','CG1514','CG15247','CG15293','CG15313','CG15316','CG15316-like', \
'CG15323','CG15327','CG15327-like','CG15336','CG15527','CG1561','CG1561-like','CG15717','CG15793-like', \
'CG15890','CG1594-like','CG1636','CG1657','CG1673','CG16743','CG16752','CG1677','CG16772','CG1689-like', \
'CG16926','CG16985','CG17012','CG1703','CG17075','CG17108','CG17234','CG17239','CG17255-like','CG17361', \
'CG1737-like','CG17376','CG17404','CG1749','CG1771-like','CG1780-like','CG1799-like','CG18067','CG18125', \
'CG1817-like','CG18265','CG18266','CG18341','CG18358','CG1839-like','CG18418','CG1847-like','CG18543', \
'CG1910','CG1950','CG1961','CG1979','CG2056','CG2059','CG2079-like','CG2145','CG2145-like','CG2150', \
'CG2217','CG2262-like','CG2446','CG2556','CG2574','CG2577','CG2662','CG2662-like','CG2829','CG2854', \
'CG2909','CG2930','CG2947','CG2984-like','CG2995','CG3004','CG3021','CG3038','CG30497','CG3071','CG3077', \
'CG3085','CG31217','CG3126-like','CG31686','CG31751','CG31826','CG31909','CG32409','CG32541','CG32549', \
'CG32582','CG32611','CG32613-like','CG32626','CG32629','CG32635','CG32637','CG32647','CG32654','CG32654-like', \
'CG32666','CG32683','CG32683-like','CG32688-like','CG32690','CG32694','CG32699','CG32700','CG32709-like', \
'CG32712','CG32732-like','CG32779','CG32789','CG32790','CG33080','CG33080-like','CG34104','CG34376','CG3476', \
'CG3483','CG3509','CG3568','CG3568-like','CG3585','CG3588','CG3592','CG3595-like','CG3600','CG3603','CG3630', \
'CG3652','CG3665-like','CG3683','CG3690','CG3699','CG3703','CG3704','CG3706','CG3708','CG3713','CG3726', \
'CG3777','CG3831','CG3918-like','CG3973','CG3975','CG4095','CG4101','CG4420','CG4420-like','CG4523-like', \
'CG4532','CG4570','CG4593','CG4593-like','CG4607','CG4645','CG4699','CG4757','CG4766-like','CG4865','CG4973', \
'CG5014-like','CG5045','CG5171','CG5177','CG5222','CG5254','CG5261','CG5273','CG5276','CG5321','CG5334', \
'CG5446','CG5539','CG5541','CG5565','CG5599','CG5662','CG5669','CG5757','CG5765','CG5773','CG5919','CG5937', \
'CG5976','CG5989','CG6036','CG6094','CG6121-like','CG6130','CG6133','CG6170-like','CG6227','CG6255','CG6332', \
'CG6459','CG6467','CG6525','CG6659','CG6687','CG6789','CG6971','CG6978','CG6980','CG6981','CG6999','CG7093', \
'CG7107-like','CG7219','CG7251','CG7387','CG7484','CG7840','CG7860','CG7920','CG7953','CG7966','CG8310', \
'CG8326','CG8334','CG8557','CG8564','CG8675','CG8909','CG8949','CG8952','CG9080','CG9114','CG9123','CG9125', \
'CG9135','CG9164','CG9245','CG9273','CG9283','CG9314','CG9355-like','CG9437','CG9518','CG9531','CG9533-like', \
'CG9571','CG9617','CG9631','CG9649','CG9650','CG9650-like','CG9723','CG9822','CG9897','CG9904','CG9907-like', \
'CG9915','CG9919','CG9928','chb','Chd1','Chi','CHIP','CHLD3','chm','Cht2','Cht3','Cht4','ci','cib','cid','cin', \
'CkIIbeta','Clic','cn','cnk','cona','cos','CoVa','COX4','COX4-2','COX5A','COX5B-1','COX5B-2','COX6A-1', \
'COX6A-2','COX6A-3','COX6B','COX6C','COX7A','COX7C','COX8','Cp15','Cp16','Cp18','Cp36','Cpr11B','Cpr62Ba', \
'Cpr76Ba','Crag','Crg-1','Crk','crm','croc','crq','Crtp','Cry','Crz','Csp','csw','cv','CycB3','cyp','Cyp1', \
'Cyp12a4','Cyp12a5','Cyp18a1','Cyp28a1','Cyp4d1','Cyp4D10','Cyp4d2','Cyp4g15','Cyp4s3','Cyp6a2','Cyp6a23', \
'Cyp6a8','Cyp6g1','Cyp9f2','cype','D','da','dah','dally','dbe','Dbp45A','Dbp73D','dc','Dcr-1','Ddc','Ddx1', \
'DebB','dec-1','deltaTry','desat 2','desatF','Dfd','Dhc-Yh3','disco','dj','Dlc90F','dlg1','dm','DnaJ-1', \
'DNApol-alpha180','dnc','Dnz1','Doa','dod','Dot','Dox-A2','dp','dpp','Dpr8','dpr8','Dpt','DptB','Dr','Dref', \
'drk','drl','dro1','dro2','dro3','dro4','dro5','dro6','Drs','dsh','Dsk','Dsor1','Dsp1','dy','eag','Edg78E', \
'Edg84A','Edg91','EG:30B8.6','Egfr','egh','Eh','eIF-4E','eIF4E-5','eIF4E-7','Eig71Ea','Eig71Eb','Eig71Ed', \
'Eig71Ee','Eig71Eg','Eig71Eh','Eig71Ei','Eig71Ej','Eig71Ek','Eip55E','elav','emb','emc','en','Eno','Ephrin', \
'esc','esg','Est-5A','Est-5B','Est-5C','Est-6','Est-A','Est-P','esterase 6c','esterase 7', \
'etaTry','Ets97D','eve','ewg','Exp6','exu1','exu2','ey','eya','f','faf','fan','Fas2','Fbp2','fd59A','fd64A', \
'fd96Ca','fd96Cb','Fdxh','Fer3','Fer3HCH','fh','fidipidine','FK506-bp2','fliI','flw','fne','fru','Fsh','Fst','ftz', \
'fu','futsch','fw','fy','fz4','fzo','g','G6pd','GA10505','GA18781','gammaCop','gammaTry','gammaTub23C','gammaTub37C', \
'Gapdh1','Gapdh2','gatA','gbb','gce','gcm','gd','gdl-ORF39','Gel','Gl','GlcAT-I','Gld','glu','Glut3','GNBP1','GNBP2', \
'GNBP3','gnu','Gpdh','Gpi1','gprs','Gr2a','Gr47A','Gr63a','grau','grk','GstD1','gt','gw','Gyc32E','gypsypol','h','H1', \
'H2BFQ','hb','Hcf','HDAC4','HDAC6','her','HERC2','Hex','Hex-A','Hex-C','HexA','hgo','His2A','His2Av','His3','His4', \
'histone H3.3','hk','Hlc','HLHm7','HLHmbeta','HLHmdelta','HmgD','Hmr','Hmu','Hn','Hop','hop','Hr4','Hrb87F','Hsc70-1', \
'Hsc70-3','Hsc70-4','Hsc70Cb','Hsf','Hsp23','Hsp26','Hsp27','Hsp68','Hsp70Bb','Hsp83','HSPA1L','htl','hug','hyd', \
'hydra','I-2','I-t','Idgf1','Idgf3','Idgf4','ifc','IM23','Imp','ImpE2','inaD','ind','InR','inx6','iotaTry','Ip259', \
'Irbp','Iris','Irp-1B','Ixl','Jafrac1','janA','janB','jgw','Jhe2','Jheh3','Jon99Cii','Jon99Ciii','jumu','KaiRIA', \
'kappaTry','kek1','kek2','kin17','kl-2','kl-5','kl2','Klp3A','Klp54D','Klp61F','klp61F','Klp64D','Klp67A','Klp68D', \
'Klp98A','kni','KP78a','KP78b','Kr','kraken','ksr','Ku80','kuz','kz','lab','Lam','lambdaTry','LanB2','Lar','lbe','Lcp1', \
'Lcp2','Lcp3','Lcp4','lcs','Lgp1','Lgp3','Lhr','lic','Lim1','Lim3','lin','Lsp1beta','Lsp1gamma','LysD','LysP','LysX', \
'lz','m1','m2','m6','mago','mal','malpha','mam','Marf','mas','Mct1','Mdr50','Mdr65','MED22','MED24','MED26','Mef2', \
'mei-218','mei-41','mei-9','mei-P22','mei-P26','mei-S332','mei-W68','Met75Ca','Met75Cb','mew','mex1','Mgstl-Psi', \
'Mhc','Mi-2','mira','Mlc1','mle','Mlh1','mmd','mof','moj30','moj9','mre11','mRpL12','mRpL14','msl-1','msl-2','msl-3', \
'msps','Mst26aa','Mst26ab','Mst57Da','Mst57Db','Mst57Dc','Mst84Da','Mst84Db','Mst84Dc','Mst84Dd','Mst85C','Mst87F', \
'Mst98Ca','Mst98Cb','mtacp1','mth','mthl8','Mtk','MtnA','mtnA','MtnB','mtrm','mud','mus101','mus205','mus304','mus308', \
'mus81','mys','N','nAcRalpha-7E','ncd','ND75','nej','Nelf-E','NetA','NetB','neur','NFAT','ng1','ng2','ng3','ng4','NHP2', \
'ninaA','ninaE','NiPp1','Nipsnap','nmd','nmdyn-D6','Nmt','nod','noi','nompA','nonA','nop5','Nop56','nos','NP15.6','npf', \
'Nrg','Nrk','Ntf-2r','NUCB1','nullo','Numb','Nup154','Nup160','Nurf-38','O-fut2','Obp56a','Obp56b','Obp56c','Obp56d', \
'Obp56e','Obp56f','Obp56g','Obp56h','Obp56i','Obp57d','Obp99a','Obp99c','Obp99d','ocn','Octbeta2R','Odc1','Odc2','okr', \
'ome','or22','or22a','or22b','Or33B','Or33b','ord','org-1-5f','Ork1','os','Os-E','Os-F','osa','Osbp','osk','Ost48','otu', \
'ovo','P5cr','pa','Pabp2','pan','para','pb','Pc','Pcl','Pcmt','pcx','Pdsw','Pen','per','Pgd','Pgi','Pgk','Pglym78','Pgm', \
'PGRP-SA','PGRP-SB1','PGRP-SB2','PGRP-SC1a','PGRP-SC1b','PGRP-SC2','PGRP-SD','ph-d','ph-p','PhKgamma','pho','Pi3K92E','pim', \
'pit','piwi','Pka-R1','Pkcdelta','Pkd2','ple','plexA','plu','Pms2','pn','pnr','poe','pol-like','polo','Pp1-13C','Pp1-87B', \
'Pp1-Y1','Pp2B-14D','PpD5','PpN58A','PPP4R2r','Prat','Prm','Pros26','Pros28.1','Pros28.1A','Pros28.1B','Pros35','Prosalpha3T', \
'Prosalpha6','Prosalpha6T','Prosbeta2','Prosbeta3','Prp18','pseudo-tra','Psf3','psh','psiI','psiII','PSR','Ptp10D','Pu', \
'px','qm','R','r','r-l','rab3-GEF','Rab7','rad','Rad23','rad50','RanGap','Ras64B','Ras85D','Rbf2','Rbp4','rdgA','rec', \
'RecQ5','Rel','Rep4','Ret','RfaBp','rg','Rh4','rha','rhi','RhoGAP15B','rl','Rlb1','Rlc1','Ro-2','robl','rok','Ror','RPA2', \
'Rpd3','Rph','RpI1','RpI135','RpII140','RpII15','RpIII128','RpL14','RpL22','RpL23','RpL32','RpL40','RpS15Ab','RpS27A','RpS3', \
'Rrp1','run','rut','rux','ry','sala','Sara','sax','sbr','sc','Scamp','scpr-C','Scr','Sd-RanGAP','SdhB','sec10','sens','Sep1', \
'Ser','Ser12','Ser7','sesB','Set','sev','sfl','sgl','Sgs1','Sgs3','Sgs4','Sgs5','Sgs7','Sgs8','Sh','sh','shark','Shc','shi', \
'sima','sina','SIP3','siren1','sisA','skl','sktl','Slip1','slp1','slp2','sls','sm','smc1','Smr','sna','snf','snf1A','SOD', \
'Sod','Sod2','Sodh-1','sog','sop','Sos','Sox102F','Sp1','Spase25','spen','spn-A','spn-B','spn-D','spri','Sptr','Spx','sqh', \
'Sr-CI','Sr-CIII','Srp19','Sry-alpha','Sry-beta','Sry-delta','Ssl1','Ssrp','sta','stan','Stellate 2','SteXh','stil','sub', \
'sunz','SuUR','swm','sws','Sxl','Syb','Syx4','T-cp1','T3dh','Taf8','Takr86C','Tango4','tara','Tbh','Tbp','Tbp-1','tefu', \
'Tehao','Ten-a','Ten-m','TepI','TepII','TepIV','Tes100','Tes101','Tes104','Tes105','Tes106','Tes107','Tes109','Tes110', \
'Tes112','Tes113','Tes114','Tes115','Tes118','Tes134','Tes14','Tes154','Tes33','TfIIA-S','TfIIB','TfIIS','TH1','Thd1', \
'thetaTry','Thor','thr','tim','tin','Tis11','tko','tld','tll','tok','Toll-7','Tollo','Tom','Top2','tor','torp4a','TotA', \
'toy','Tpi','Tps1','tra','tra2','trc','Tre','Tre1','Treh','Trf','Trf2','tRNA:Y1:22Fa','trol','tropomyosin isoform 9A', \
'trp','trpl','trr','Ts','tsg','Tsp33B','Tsp3A','tsr','ttk','tud','Uba1','UBC9','UbcD4','Ubx','UGP','unc-119','unc-13', \
'Unc-76','UP','up','Updo','Uro','v','vas','vg','Vha68-1','Vinc','vir','vlc','vls','Vm26Aa','Vm26Ab','Vm32E','Vm34Ca', \
'w','wapl','waw','wbl','wdn','wds','wg','wgn','wit','woc','X11Lbeta','Xdh','xdh','Xpac','Xpd','y','Yp1','Yp2','Yp3', \
'z','Z600','ZAMenv','Zeelin1','zen','zen2','zetaTry','zfh2','Zw']

        '''automaitic research pipe00_00_Search_Species_x_Genes.py '''
        self.drosophilaGoodGenes2 = ['Adh', 'Adh-1', 'Adh-2', 'Adh1', 'Adh2', 'AdhPsi', 'Adhr', 'AdhT', \
                       'amyRel', 'amy','Amy-d','Amy-p','amy-d','Amy1','Amy2', \
                      'bcd', 'cdc6', 'Cp36', 'Cp16', \
                      'CG33080','CG12772', 'CG15327', 'CG11802', 'CG10996', 'CG3568', \
                      'Ddc', 'Dhc-Yh3', 'Eno', 'ey', 'ftz', 'G6pd', 'Gpdh', 'Gapdh2', \
                      'Hex-A', 'Hex-C', 'Hmr', 'hydra',\
                       'N', 'sqh', \
                       'mth',  'Os-E', 'Os-F', 'Pgd','Pgi',  \
                      'Mlc1', 'Acp53Ea', 'w', 'v', 'Ku80', 'Acp26Aa', 'Sxl', 'per', 'run', 'rux', 'Tpi','tra', 'y']

        self.drosophilaGoodGenes = ['Acp26Aa', 'Acp53Ea', \
                       'Adh', 'Adh-1', 'Adh-2', 'Adh1', 'Adh2', 'AdhPsi', 'Adhr', 'AdhT', \
                       'bcd', 'CG10996', 'CG11802', 'CG12772', 'CG15327', 'CG33080', 'CG3568', \
                       'cdc6', 'Cp16', 'Cp36', 'Ddc', 'Dhc-Yh3', 'Eno', 'ey', 'ftz', 'G6pd', 'Gapdh2', 'Gpdh', \
                       'Hex-A', 'Hex-C', 'Hmr', 'hydra', 'Ku80', 'Mlc1', 'mth', 'N', 'Os-E', 'Os-F', \
                       'per', 'Pgd', 'Pgi', 'run', 'rux', 'sqh', 'Sxl', 'Tpi', 'tra', 'v', 'w', 'y']


        self.drosophilaSomeGoodGenes = ['Acp26Aa', 'Acp53Ea', \
                       'bcd', 'CG10996', 'CG11802', 'CG12772', 'CG15327', 'CG33080', 'CG3568', \
                       'cdc6', 'Cp16', 'Cp36', 'Dhc-Yh3', 'Eno', 'ey', 'ftz', 'G6pd', 'Gapdh2', \
                       'Hex-A', 'Hex-C', 'Hmr', 'hydra', 'Ku80', 'Mlc1', 'mth', 'N', 'Os-E', 'Os-F', \
                       'Pgd', 'run', 'rux', 'sqh', 'Sxl', 'Tpi', 'tra', 'v', 'w', 'y']



    def databaseInfo(self):
        
        handle = Entrez.einfo()
        record = Entrez.read(handle)
        
        return record["DbList"]


        '''
        ['pubmed', 'protein', 'nucleotide', 'nuccore', 'nucgss', 'nucest',
         'structure', 'genome', 'books', 'cancerchromosomes', 'cdd', 'gap',
         'domains', 'gene', 'genomeprj', 'gensat', 'geo', 'gds', 'homologene',
         'journals', 'mesh', 'ncbisearch', 'nlmcatalog', 'omia', 'omim', 'pmc',
         'popset', 'probe', 'proteinclusters', 'pcassay', 'pccompound',
         'pcsubstance', 'snp', 'taxonomy', 'toolkit', 'unigene', 'unists']
        '''
   
    def adhCorrection(self, stri):
        stri = stri.replace('pseudoobscura subsp bogotana', 'pseudoobscura bogotana')
        return 
    
    def mnemonic(self, species):
        species = species.replace("-","_")
        
        if species == 'albomicans':
            return 'albo'
        elif species == 'americana':
            return 'am'
        elif species == 'americana_americana'  or species == 'americanaamericana':
            return 'am-am'
        elif species == 'americana_texana' or species == 'americanatexana':
            return 'am-tx'
        elif species == 'ananassae':
            return 'ana'
        elif species == 'angor':
            return 'angor'
        elif species == 'antillea':
            return 'antil'
        elif species == 'arawakana':
            return 'araw'
        elif species == 'arawakana_arawakana':
            return 'araw2'
        elif species == 'arawakana_kittensis':
            return 'arawkit'
        elif species == 'arizonae':
            return 'ariz'
        elif species == 'bipectinata':
            return 'bipect'
        elif species == 'caribiana':
            return 'carib'
        elif species == 'grimshawi':
            return 'grimx'
        elif species == 'innubila':
            return 'innub'
        elif species == 'kikkawai':
            return 'kik'
        elif species == 'malerkotliana':
            return 'malerk'
        elif species == 'malerkotliana_malerkotliana':
            return 'malerk2'
        elif species == 'mayaguana':
            return 'maya'
        elif species == 'mauritiana':
            return 'mauri'
        elif species == 'melanogaster':
            return 'mela'
        elif species == 'nigrodunni':
            return 'nigrod'
        elif species == 'miranda':
            return 'mira'
        elif species == 'mojavensis':
            return 'moja'
        elif species == 'nasuta':
            return 'nasut'
        elif species == 'nigrodunni':
            return 'nigro'
        elif species == 'novamexicana':
            return 'novmex'
        elif species == 'parisiena':
            return 'paris'
    
        elif species == 'paulistorum':
            return 'paul'
        elif species == 'persimilis':
            return 'persi'
        elif species == 'pseudoobscura':
            return 'pseud'
        elif species == 'pseudoobscura_bogotana':
            return 'pbogo'
        elif species == 'pseudoobscura_pseudoobscura':
            return 'pseud2'
        elif species == 'recens':
            return 'recens'
        elif species == 'santomea':
            return 'stom'
        elif species == 'sechellia':
            return 'sec'
        elif species == 'similis':
            return 'simi'
        elif species == 'similis_similis':
            return 'simi2'
        elif species == 'simulans':
            return 'sim'
        elif species == 'subobscura':
            return 'subo'
        elif species == 'sturtevanti':
            return 'sturt'
        elif species == 'straubae':
            return 'straub'
        elif species == 'teissieri':
            return 'teis'
        elif species == 'triauraria':
            return 'tria'
        elif species == 'virilis':
            return 'viri'
        elif species == 'willistoni':
            return 'wil'
        elif species == 'yakuba':
            return 'yak'
        
        return species
        
    def labels(self, speciesList):
        lista  = []
        for species in speciesList:
            species = species.replace("-","_")
            
            if species == 'albomicans':
                lista.append('albo')
            elif species == 'americana':
                lista.append('am')
            elif species == 'americana_americana':
                lista.append('am-am')
            elif species == 'americana_texana':
                lista.append('am-tx')
            elif species == 'ananassae':
                lista.append('ana')
            elif species == 'angor':
                lista.append('angor')
            elif species == 'antillea':
                lista.append('antil')
            elif species == 'arawakana':
                lista.append('araw')
            elif species == 'arawakana_arawakana':
                lista.append('araw2')
            elif species == 'arawakana_kittensis':
                lista.append('arawkit')
              
            elif species == 'arizonae':
                lista.append('ariz')
            elif species == 'bipectinata':
                lista.append('bipect')
            elif species == 'bipectinata':
                lista.append('bipect')
            elif species == 'caribiana':
                lista.append('carib')
            elif species == 'grimshawi':
                lista.append('grimx')
            elif species == 'innubila':
                lista.append('innub')
            elif species == 'mauritiana':
                lista.append('mauri')
            elif species == 'mayaguana':
                lista.append('maya')

            elif species == 'malerkotliana':
                lista.append('malerk')
            elif species == 'malerkotliana_malerkotliana':
                lista.append('malerk2')

            elif species == 'melanogaster':
                lista.append('mela')
            elif species == 'miranda':
                lista.append('mira')
            elif species == 'mojavensis':
                lista.append('moja')
            elif species == 'nigrodunni':
                lista.append('nigrod')
            elif species == 'novamexicana':
                lista.append('novmex')
            elif species == 'parisiena':
                lista.append('paris')
            elif species == 'persimilis':
                lista.append('persi')
            elif species == 'pseudoobscura':
                lista.append('pseud')
            elif species == 'pseudoobscura_bogotana':
                lista.append('pbogo')
            elif species == 'pseudoobscura_pseudoobscura':
                lista.append('pseud2')
            elif species == 'recens':
                lista.append('recens')
            elif species == 'santomea':
                lista.append('stom')
            elif species == 'similis':
                lista.append('simi')
            elif species == 'similis_similis':
                lista.append('simi2')
            elif species == 'sechellia':
                lista.append('sec')
            elif species == 'simulans':
                lista.append('sim')
            elif species == 'sturtevanti':
                lista.append('sturt')
            elif species == 'straubae':
                lista.append('straub')
                
            elif species == 'subobscura':
                lista.append('subo')
            elif species == 'sturtevanti':
                lista.append('sturt')             
            elif species == 'teissieri':
                lista.append('teis')
            elif species == 'triauraria':
                lista.append('tria')                
            elif species == 'virilis':
                lista.append('viri')
            elif species == 'willistoni':
                lista.append('wil')
            elif species == 'yakuba':
                lista.append('yak')
            else:
                lista.append(species)   
                
        return lista


