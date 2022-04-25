#Initally start with eQTL to work out how to do it. 3 steps required:
#* get all of the probe names in the file
#* create a folder with the probe name
#* Create list of SNPs associated with the probe
#* In each folder save the individual level data for phenotypes and genottypes for the SNPs

#combine the individual level data down to one file so can be used to find each data for probe
cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data
#combine all of the individual level data for each chromosome without header
awk -F, FNR!=1 chr* > all.txt
#add the header to the data
( echo -e "CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8	sample9	sample10	sample11	sample12	sample13	sample14	sample15	sample16	sample17	sample18	sample19	sample20	sample21	sample22	sample23	sample24	sample25	sample26	sample27	sample28	sample29	sample30	sample31	sample32	sample33	sample34	sample35	sample36	sample37	sample38	sample39	sample40	sample41	sample42	sample43	sample44	sample45	sample46	sample47	sample48	sample49	sample50	sample51	sample52	sample53	sample54	sample55	sample56	sample57	sample58	sample59	sample60	sample61	sample62	sample63	sample64	sample65	sample66	sample67	sample68	sample69	sample70	sample71	sample72	sample73	sample74	sample75	sample76	sample77	sample78	sample79	sample80	sample81	sample82	sample83	sample84	sample85	sample86	sample87	sample88	sample89	sample90	sample91	sample92	sample93	sample94	sample95	sample96	sample97	sample98	sample99	sample100	sample101	sample102	sample103	sample104	sample105	sample106	sample107	sample108	sample109	sample110	sample111	sample112	sample113	sample114	sample115	sample116	sample117	sample118	sample119	sample120	sample121	sample122	sample123	sample124	sample125	sample126	sample127	sample128	sample129	sample130	sample131	sample132	sample133	sample134	sample135	sample136	sample137	sample138	sample139	sample140	sample141	sample142	sample143	sample144	sample145	sample146	sample147	sample148	sample149	sample150	sample151	sample152	sample153	sample154	sample155	sample156	sample157	sample158	sample159	sample160	sample161	sample162	sample163	sample164	sample165	sample166	sample167	sample168	sample169	sample170	sample171	sample172	sample173	sample174	sample175	sample176	sample177	sample178	sample179	sample180	sample181	sample182	sample183	sample184	sample185	sample186	sample187	sample188	sample189	sample190	sample191	sample192	sample193	sample194	sample195	sample196	sample197	sample198	sample199	sample200	sample201	sample202	sample203	sample204	sample205	sample206	sample207	sample208	sample209	sample210	sample211	sample212	sample213	sample214	sample215	sample216	sample217	sample218	sample219	sample220	sample221	sample222	sample223	sample224	sample225	sample226	sample227	sample228	sample229	sample230	sample231	sample232	sample233	sample234	sample235	sample236	sample237	sample238	sample239	sample240	sample241	sample242	sample243	sample244	sample245	sample246	sample247	sample248	sample249	sample250	sample251	sample252	sample253	sample254	sample255	sample256	sample257	sample258	sample259	sample260	sample261	sample262	sample263	sample264	sample265	sample266	sample267	sample268	sample269	sample270	sample271	sample272	sample273	sample274	sample275	sample276	sample277	sample278	sample279	sample280	sample281	sample282	sample283	sample284	sample285	sample286	sample287	sample288	sample289	sample290	sample291	sample292	sample293	sample294	sample295	sample296	sample297	sample298	sample299	sample300	sample301	sample302	sample303	sample304	sample305	sample306	sample307	sample308	sample309	sample310	sample311	sample312	sample313	sample314	sample315	sample316	sample317	sample318	sample319	sample320	sample321	sample322	sample323	sample324	sample325	sample326	sample327	sample328	sample329	sample330	sample331	sample332	sample333	sample334	sample335	sample336	sample337	sample338	sample339	sample340	sample341	sample342	sample343	sample344	sample345	sample346	sample347	sample348	sample349	sample350	sample351	sample352	sample353	sample354	sample355	sample356	sample357	sample358	sample359	sample360	sample361	sample362	sample363	sample364	sample365	sample366	sample367	sample368	sample369	sample370	sample371	sample372	sample373	sample374	sample375	sample376	sample377	sample378	sample379	sample380	sample381	sample382	sample383	sample384	sample385	sample386	sample387	sample388	sample389	sample390	sample391	sample392	sample393	sample394	sample395	sample396	sample397	sample398	sample399	sample400	sample401	sample402	sample403	sample404	sample405	sample406	sample407	sample408	sample409	sample410	sample411	sample412	sample413	sample414	sample415	sample416	sample417	sample418	sample419	sample420	sample421	sample422	sample423	sample424	sample425	sample426	sample427	sample428	sample429	sample430	sample431	sample432	sample433	sample434	sample435	sample436	sample437	sample438	sample439	sample440	sample441	sample442	sample443	sample444	sample445	sample446	sample447	sample448	sample449	sample450	sample451	sample452	sample453	sample454	sample455	sample456	sample457	sample458	sample459	sample460	sample461	sample462	sample463	sample464	sample465	sample466	sample467	sample468	sample469	sample470	sample471	sample472	sample473	sample474	sample475	sample476	sample477	sample478"; cat all.txt ) > combined_individual_eQTL.txt
rm all.txt


#seperate and get summary and individual level data for each SNP associated with each probe
cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping
mkdir single_probes
cd single_probes
pwd #/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes

#get the probe names
cut -d',' -f1 /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/map_to_rsId/sorted_eqtl_results_p0.01.csv | sort | uniq
#prints list of the unique probes in the file

#copy all of the single probes to the scratch to save space and speed up this and later steps.
cp -R /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/clumping/single_probes ./
cp /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/map_to_rsId/sorted_eqtl_results_p0.01.csv ./
cp /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/combined_individual_eQTL.txt ./

#Create folder for each probe and Create list of SNPs associated with the probe
#Run on server
#!/bin/bash --login
#$ -cwd               # Run job from current directory
for OUTPUT in $(cut -d',' -f1 /mnt/iusers01/bk01/m20349jk/scratch/single_probes/sorted_eqtl_results_p0.01.csv | sort | uniq)
do
  mkdir ${OUTPUT}
  cd ${OUTPUT}
  awk -F , -v myvar=$OUTPUT '$1 == myvar ' /mnt/iusers01/bk01/m20349jk/scratch/single_probes/sorted_eqtl_results_p0.01.csv > associated_SNPs.csv
  awk -F"," '{print $5} ' associated_SNPs.csv > SNPs.csv #create file of SNPs
  fgrep -w -f SNPs.csv /mnt/iusers01/bk01/m20349jk/scratch/single_probes/combined_individual_eQTL.txt > individual_SNPs.csv #get list of indiv SNPs data for the SNPs
  sed -i "1s/^/$(head -n1 /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/combined_individual_eQTL.txt)\n/" individual_SNPs.csv #add heading to the SNP table
  echo "${OUTPUT} completed"
  cd ..
done



qsub sep_by_probe.txt
rm -r GeneID

