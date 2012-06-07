#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <algorithm>
#include <time.h>
 

struct Frame
{
    int ObjID1;
    int SObjID;
    int NodeID1;
    //int NormObjID;
               
};


int struct_cmp(const void *a, const void *b)
{
    struct Frame *ia = (struct Frame *)a;
    struct Frame *ib = (struct Frame *)b;
    return (int) (ia->NodeID1 - ib->NodeID1);
}

int struct_cmpObj(const void *x, const void *y)
{
    struct Frame *ix = (struct Frame *)x;
    struct Frame *iy = (struct Frame *)y;
    return (int) (ix->SObjID - iy->SObjID);
}



using namespace std;

int main()
{
    int filecounter1, filecounter2;
    char st1[32] = "halos_130_subhalo_";
    char st2[32] = "halos_138_subhalo_";
    char ext[10] = ".cosmo";
    char filename [50] = {'\0'};
    float pos1, vel1, pos2, vel2, pos3, vel3, tag1; //Positions, velocities and subhalo_tags
    int tag2; //Particle tag
    int high_sub1 = 0, high_sub2 = 0;
    int increment = 0;
    int Halo_count1 = 0, Halo_count2 = 0, Global_halo;
    int count=0, countN=0, countO=0, Cnode, Cobj, count1=0, countN2=0, countO2=0, Cnode1, Cobj1, wait;
    clock_t start1, start2;
    double duration;
    vector<int> Halo1;
    vector<int> Halo2;
    vector<int> VSobj;
    vector<int> VObj;
    vector<int> VNode;


   

  // Reads binary files for first time step.
  for (filecounter1 = 1; filecounter1 < 600; filecounter1++)
    {
        sprintf_s(filename, "%s%d%s", st1, filecounter1, ext);

        ifstream file (filename, ios::in|ios::binary);
        int iter = 0; //To categorise the particle parameters according to their position in the file

        if(file.is_open())
        {
       
            while(!file.eof())
            {

                if ((iter % 8) == 0)
                {
                    file.read((char *)(&pos1), sizeof(pos1));
                }
       

                else if ((iter % 8) == 1)
                {
                    file.read((char *)(&vel1), sizeof(vel1));
                }
   

                else if ((iter % 8) == 2)
                {
                    file.read((char *)(&pos2), sizeof(pos2));
                }


                else if ((iter % 8) == 3)
                {
                    file.read((char *)(&vel2), sizeof(vel2));
                }


                else if ((iter % 8) == 4)
                {
                    file.read((char *)(&pos3), sizeof(pos3));
                }


                else if ((iter % 8) == 5)
                {
                    file.read((char *)(&vel3), sizeof(vel3));
                }


                else if ((iter % 8) == 6) //Subhalo_tag
                {
                    file.read((char *)(&tag1), sizeof(tag1));
                                                            //cout << "subhalo binary " << tag1 << endl;
                    tag1 = tag1 + increment; // Assigns a global subhalo_tag
                    if(tag1 > high_sub1) // Records the highest subhalo_tag in this file
                    {
                        high_sub1 = tag1;
                    }
                    VSobj.push_back(tag1);
                                            //Ftime1[countO].SObjID = tag1;
                    VObj.push_back(filecounter1); //Assign halo ID
                                            //Ftime1[countO].ObjID1 = filecounter;
                                            //cout << "Subhalo IDs " << VSobj[countO] << "Halo ID " << VObj[countO] << endl;
                    countO++;
                }
       

                else if ((iter % 8) == 7) //Particle_tag
                {
                    file.read((char *)(&tag2), sizeof(tag2));
                    VNode.push_back(tag2);
                                            //Ftime1[countN].NodeID1 = tag2;
                                            //cout << "Particle ID " << VNode[countN] << endl;
                    countN++;
                }

                iter++;

            }

            Halo_count1++; // Count the number of Halos in time i.
        }
        file.close();
        file.clear();
           
        increment = (high_sub1 + 1); //Assign value of higest subhalo_tag to continue numbering
    }
    Global_halo = Halo_count1;
    cout<<"Global Halo "<<Global_halo<<endl;
   
    Frame* Ftime1;
      Ftime1 = new Frame[countO];
    for(int sub = 0; sub < countO; sub++) //Trasfer info to vectors
    {
        Ftime1[sub].SObjID = VSobj[sub];
        Ftime1[sub].ObjID1 = VObj[sub];
        Ftime1[sub].NodeID1 = VNode[sub];
                                            //cout << "Subhalo " << Ftime1[sub].SObjID << "Halo " << Ftime1[sub].ObjID1 << "Node " << Ftime1[sub].NodeID1 << endl;
    }
    VSobj.clear();
    VObj.clear();
    VNode.clear();

    cout<<"Ftime1 Counts"<<countO<<" "<<countN<<endl;
    cout<<"Sizeof"<<(sizeof(Ftime1)/sizeof(Ftime1[0]))<<endl;
    cout << "Highsub " << high_sub1 << endl;

 
    vector<int> HaloVol1(high_sub1 + 1);

    size_t Ftime1_len = sizeof(Ftime1) / sizeof(struct Frame);


    start1 = clock();
    qsort(Ftime1,countO, 12, struct_cmpObj); //Ftime1_len = countO & sizeof(struct Frame) = total bytes of all members
    duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
    //cout << "Duration1 : " << duration << endl;



   for(int Ncount = 0; Ncount < countO; Ncount++)
   {
       Halo1.push_back(Ftime1[Ncount].SObjID);      
       HaloVol1[Ftime1[Ncount].SObjID]++; //For counting number of particles in halo
   }


   start2 = clock();
   qsort(Ftime1,countO, 12, struct_cmp);
   duration = (clock() - start2); // (double) CLOCKS_PER_SEC;
   //cout << "Duration2 : " << duration << endl;
   cout << "Sorted again!" << endl;


  


   // Reads binary files for second time step.
  increment = 0;

  for (filecounter2 = 1; filecounter2 < 600; filecounter2++)
    {
        sprintf_s(filename, "%s%d%s", st2, filecounter2, ext);

        ifstream file2 (filename, ios::in|ios::binary);
        int iter2 = 0; //To categorise the parameters according to their position in the file

        if(file2.is_open())
        {
       
            while(!file2.eof())
            {

                if ((iter2 % 8) == 0)
                {
                    file2.read((char *)(&pos1), sizeof(pos1));
                }
       

                else if ((iter2 % 8) == 1)
                {
                    file2.read((char *)(&vel1), sizeof(vel1));
                }
   

                else if ((iter2 % 8) == 2)
                {
                    file2.read((char *)(&pos2), sizeof(pos2));
                }


                else if ((iter2 % 8) == 3)
                {
                    file2.read((char *)(&vel2), sizeof(vel2));
                }


                else if ((iter2 % 8) == 4)
                {
                    file2.read((char *)(&pos3), sizeof(pos3));
                }


                else if ((iter2 % 8) == 5)
                {
                    file2.read((char *)(&vel3), sizeof(vel3));
                }


                else if ((iter2 % 8) == 6)
                {
                    file2.read((char *)(&tag1), sizeof(tag1));
                                                            //cout << "subhalo binary " << tag1 << endl;
                    tag1 = tag1 + increment; // Assigns a global subhalo_tag
                    if(tag1 > high_sub2) // Records the highest subhalo_tag in this file
                    {
                        high_sub2 = tag1;
                    }
                    VSobj.push_back(tag1);
                                                            //Ftime2[countO2].SObjID = tag1;
                    VObj.push_back(filecounter2);
                                                            //Ftime2[countO2].ObjID1 = filecounter;
                                                            //cout << "Subhalo IDs " << Ftime2[countO2].SObjID << endl;
                    countO2++;
                }
       

                else if ((iter2 % 8) == 7)
                {
                    file2.read((char *)(&tag2), sizeof(tag2));
                    VNode.push_back(tag2);
                    countN2++;
                }

                iter2++;

            }

            Halo_count2++; // Count the number of halos in time i+1.
        }
        file2.close();
        file2.clear();
       
   
        increment = (high_sub2 + 1); //Assign value of highest subhao_tag to continue numbering
    }
 


    Frame* Ftime2;
    Ftime2 = new Frame[countO2];
    for(int sub2 = 0; sub2 < countO2; sub2++)
    {
        Ftime2[sub2].SObjID = VSobj[sub2];
        Ftime2[sub2].ObjID1 = VObj[sub2];
        Ftime2[sub2].NodeID1 = VNode[sub2];
    }
    VSobj.clear();
    VObj.clear();
    VNode.clear();


   cout<<"Ftime2 Counts"<<countO2<<" "<<countN2<<endl;
   cout<<"Sizeof"<<(sizeof(Ftime2)/sizeof(Ftime2[0]))<<endl;
   cout << "Highsub " << high_sub2 << endl;


   vector<int> HaloVol2(high_sub2 + 1);

   size_t Ftime2_len = sizeof(Ftime2) / sizeof(struct Frame);

   start1 = clock();
   qsort(Ftime2,countO2, 12, struct_cmpObj);
   duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
   //cout << "Duration1 : " << duration << endl;

                                                        //Ftime2[0].NormObjID = 0;
   for(int Ncount2 = 0; Ncount2 < countO2; Ncount2++)
   {
        Halo2.push_back(Ftime2[Ncount2].SObjID);
        HaloVol2[Ftime2[Ncount2].SObjID]++;
   }

     
   start2 = clock();
   qsort(Ftime2,countO2, 12, struct_cmp);
   duration = (clock() - start2); // (double) CLOCKS_PER_SEC;
   //cout << "Duration2 : " << duration << endl;
   //cout<<"Sizeof after sort"<<(sizeof(Ftime2)/sizeof(Ftime2[0]))<<endl;
 


     int numNodes1=(sizeof(Ftime1)); ///sizeof(Ftime1[0]));
     //cout<<"Size Ftime1 Numnodes"<<numNodes1<<endl;
     int numNodes2=(sizeof(Ftime2)); ///sizeof(Ftime2[0]));
     //cout<<"Size Ftime2 Numnodes"<<numNodes2<<endl;
   
     
     vector< vector <int> > OverlapTable ((high_sub1 + 2), vector<int> (high_sub2 + 2)); // Initialise the subhalo Overlap Table
     //vector< vector <int> > OverlapTable_H ((Halo_count1 + 1), vector<int> (Halo_count2 + 1)); // Initialise the Halo Overlap Table
     register int iT (0), jT (0);


       start1 = clock();
     //Overlap detection
     while (iT<countO && jT<countO2)
     {
         /*if (iT == 933 || jT == 1551)
         {
         cout<<"In While"<<endl;
         cout<<"iT "<<iT<<"numNodes1 "<<numNodes1<<endl;
         cout<<"jT "<<jT<<"numNodes2 "<<numNodes2<<endl;
         cin >> wait;
         }
/*          if (shift_thresh > j)
          {
            temp_shift_thresh = (int)(j/2);
          }
          else
          {
            temp_shift_thresh = shift_thresh;
          }
*/


          //(t1.nodes[i].NodeID == t2.nodes[j].NodeID)   /* overlap */    //     (t1.nodes[i].NodeID == t2.nodes[j].NodeID)
          //int tempval =  (t2.nodes[j].NodeID) - (t1.nodes[i].NodeID);
            //if ( tempval <= shift_thresh)
            //if ( (tempval <= shift_thresh) && ( tempval >= 0) ) //(t1.nodes[i].NodeID == t2.nodes[j].NodeID)
            //if (t1.nodes[i].NodeID == t2.nodes[j].NodeID)  //(t1.nodes[i].NodeID == t2.nodes[j].NodeID)
      if (Ftime1[iT].NodeID1 == Ftime2[jT].NodeID1)  //( (tempval <= shift_thresh) && ( tempval >= 0) ) //(t1.nodes[i].NodeID == t2.nodes[j].NodeID)
          {
              OverlapTable[Ftime1[iT].SObjID][Ftime2[jT].SObjID]++; //Original: OverlapTable[Ftime1[iT++].ObjID1][Ftime2[jT++].ObjID1]++;
              //cout<<"An entry entered to the OVERLAPTABLE !! "<<Ftime1[iT].SObjID<<" and "<<Ftime2[jT].SObjID<<" shares node "<<Ftime1[iT].NodeID1<<endl;
              OverlapTable[high_sub1 + 1][Ftime2[jT].SObjID] = Ftime2[jT].ObjID1; // Connect subhalo to its Halo in time step i+1
              OverlapTable[Ftime1[iT].SObjID][high_sub2 + 1] = Ftime1[iT].ObjID1; // Connect subhalo to its Halo in time step i

              iT++;
              jT++;
          }

          else
         
          {

                if(Ftime1[iT].NodeID1 > Ftime2[jT].NodeID1)
                {
                    OverlapTable[high_sub1 + 1][Ftime2[jT].SObjID] = Ftime2[jT].ObjID1; // Connect subhalo to its Halo in time step i+1
                    //OverlapTable[Ftime1[iT].SObjID][high_sub2 + 1] = Ftime1[iT].ObjID1; // Connect subhalo to its Halo in time step i

                    if (jT<countO2)
                    {
                         jT++;
                    }

                }
                else
                {
                    //OverlapTable[high_sub1 + 1][Ftime2[jT].SObjID] = Ftime2[jT].ObjID1; // Connect subhalo to its Halo in time step i+1
                    OverlapTable[Ftime1[iT].SObjID][high_sub2 + 1] = Ftime1[iT].ObjID1; // Connect subhalo to its Halo in time step i

                    if (iT<countO)
                     {
                        iT++;
                     }

                }
          }
     
          /*if (iT==numNodes1 || jT==numNodes2)
          { // one frame finished, then terminate search
           iT = numNodes1;
           jT = numNodes2;
          }*/
     }
    duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
    //cout << "Duration Search : " << duration << endl;

         
   


     // Write tracking results to a file
     ofstream myfile;
     myfile.open("T1T2data_V2.txt");

     start1 = clock();
     for(int Scount1 = 0; Scount1 < (high_sub1 + 1); Scount1++)
     {
         vector <int> tempPer;
         vector <int> temp;
         vector <int> tempB;
         vector <int> tempM;
         int checkE = 0, iterB = 0, scanB = 0;

         // To check number of overlapping features in timestep i+1
         for(int Scount2 = 0; Scount2 < (high_sub2 + 1); Scount2++)
         {
             if(OverlapTable[Scount1][Scount2]>0)
             {
                tempPer.push_back(OverlapTable[Scount1][Scount2]);
                temp.push_back(Scount2);
                checkE++; // To check number of overlaps in timestep i+1
             }
                          
             else if((Scount1 == 0) && (OverlapTable[Scount1][Scount2] == 0))
             {
                 tempB.push_back(Scount2);
                 iterB++;
             }

             /*if(OverlapTable[Scount1][Scount2] > 0)
             {
                    myfile << "HaloID "<< Ftime1[Scount1].ObjID1<<" in Time199 consists of HaloID"<<Ftime2[Scount2].ObjID1<<" in Time249"<<endl;
             }*/
                                
         }


             float PerOver = 0.0;
             double var1= 0.0, var2 = 0.0;
             if((Scount1 == 0) && (iterB > 0)) //Check for BIRTH event
             {         
                 for(int birth = 0; birth < iterB; birth++)
                 {
                     scanB = 0; //Reinitialise birth check counter
                     for(int chkB = 0; chkB < (high_sub1 + 1); chkB++)
                     {
                         if(OverlapTable[chkB][tempB[birth]] > 0)
                             scanB++;
                     }
                     
                     if(scanB == 0)
                     {
                          //duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
                          //cout << "Duration Birth : " << duration << endl;
                         OverlapTable[high_sub1 + 1][tempB[birth]] = Global_halo + 1; //Use the Global_halo to record birth of a halo in timestep ti+1
                         Global_halo++; //Increment
                          myfile << "                                SubHalo ID " << tempB[birth] << " of Halo ID " << OverlapTable[high_sub1 + 1][tempB[birth]] << " is born." << endl << endl;
                     }
                 }
             }
             



             if(checkE == 0) //DEATH event
             {
                 //duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
                 //cout << "Duration Death : " << duration << endl;
                 myfile << "SubHalo ID " << Scount1 << " of Halo ID " << OverlapTable[Scount1][high_sub2 + 1] << " dies" << endl;
             }




              if(checkE>1) //SPLIT event
             {
                 //duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
                 //cout << "Duration Split : " << duration << endl;
                 myfile<<"SubHalo ID "<<Scount1<<" ("<<HaloVol1[Scount1]<<") of Halo ID " << OverlapTable[Scount1][high_sub2 + 1] << endl; //NormObj: (Scount1) else use Halo1[Scount1]
                 for(int split = 0; split < checkE; split++)
                 {
                     var1 = OverlapTable[Scount1][temp[split]];
                     var2 = HaloVol1[Scount1];
                     PerOver = (var1/var2)*100;
                     myfile<<"                     ->> ("<<PerOver<<"%) ->> SubHalo ID "<<temp[split]<<" ("<<HaloVol2[temp[split]]<<") of Halo ID " << OverlapTable[high_sub1 + 1][temp[split]] <<endl; //NormObj: temp[split] else use Halo2[Scount2]
                     cout<<"Split"<<endl;
                 }
             }


             
             if(checkE == 1) //Condition for MERGE/CONTINUATION event
             {
                 int countM =0;
                 for(int merge = 0; merge < (high_sub1 + 1); merge++)
                 {
                     if(OverlapTable[merge][temp[0]]>0)
                     {
                         tempM.push_back(merge);
                         countM++;
                     }
                 }

                 if(countM>1) //MERGE EVENT
                 {
                     //duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
                     //cout << "Duration Merge : " << duration << endl;
                     int TotOvHalo = 0, TotHalo1 = 0;
                     int Halo_IdMin = OverlapTable[tempM[0]][high_sub2 + 1];

                     for(int Cmerge = 0; Cmerge < countM; Cmerge++)
                    {
                        float PerInOver = 0.0;
                        double InVar1 = 0.0, InVar2 = 0.0;
                        InVar1 = OverlapTable[tempM[Cmerge]][temp[0]];
                        InVar2 = HaloVol1[tempM[Cmerge]];
                        PerInOver = (InVar1/InVar2)*100;
                       
                        myfile<<"SubHalo ID "<<tempM[Cmerge]<<" ("<<HaloVol1[tempM[Cmerge]]<< ") of Halo ID " << OverlapTable[tempM[Cmerge]][high_sub2 + 1] << " (" << PerInOver << "%) " <<endl; //NormObj: tempM[Cmerge] else use Halo1[tempM[Cmerge]]
                        if(OverlapTable[tempM[Cmerge]][high_sub2 + 1] < Halo_IdMin)
                        {
                            Halo_IdMin = OverlapTable[tempM[Cmerge]][high_sub2 + 1];
                        }

                        TotOvHalo = OverlapTable[tempM[Cmerge]][temp[0]] + TotOvHalo; //For percentage calculations
                        TotHalo1 = HaloVol1[tempM[Cmerge]] + TotHalo1; //For percentage calculations
                    }
                     var1 = TotOvHalo;
                     var2 = TotHalo1;
                     PerOver = (var1/var2)*100;
                     OverlapTable[high_sub1 + 1][temp[0]] = Halo_IdMin;                     
                     myfile<<"                                 }-> ("<<PerOver<<"%) }-> SubHalo ID "<<temp[0]<<" ("<<HaloVol2[temp[0]]<<") of Halo ID " << OverlapTable[high_sub1 + 1][temp[0]] <<endl; //NormObj: temp[0] else use Halo2[temp[0]]
                 }


                 else //CONTINUATION EVENT
                 {
                     //duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
                     //cout << "Duration Continuation : " << duration << endl;
                     var1 = OverlapTable[Scount1][temp[0]];
                     var2 = HaloVol1[Scount1];
                     PerOver = (var1/var2)*100; //Percentage overlap
                     OverlapTable[high_sub1 + 1][temp[0]] = OverlapTable[Scount1][high_sub2 + 1]; //Assign HaloID of timestep ti to the continuing halo in ti+1
                     myfile<<"SubHalo ID "<<Scount1<<" ("<<HaloVol1[Scount1]<<") of Halo ID " << OverlapTable[Scount1][high_sub2 + 1] << " -> ("<<PerOver<<"%) -> SubHalo ID "<<temp[0]<<" ("<<HaloVol2[temp[0]]<<") of Halo ID " << OverlapTable[high_sub1 + 1][temp[0]] <<endl; //NormObj: Scount1, temp[0]
                    //Changed Halo1[Scount1] to Scount1
                 }
             }


         myfile << endl<<endl;
         temp.clear();
     }
     duration = (clock() - start1); // (double) CLOCKS_PER_SEC;
     //cout << "Duration Classifier : " << duration << endl;


     myfile.close(); //Close file

}
