#define _XOPEN_SOURCE 500
#include <ftw.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#define fst 3 // (HaloID - 1) of the halo whose merger tree is needed. For halo 4, fst = 3. 

using namespace std;


vector<string> listFiles; // To store file names
string tempLoad;



// HaloID, SubHaloID and particleID for each particle
struct Frame
{
    int ObjID;
    int SObjID;
    int NodeID;         
};


// Velocity and position for each halo
struct Halo_param
{
    float vx;
    float vy;
    float vz;
    float px;
    float py;
    float pz;
};


// For each halo's merger info
struct Merger_info
{
    int host_node;
    int merge_node;
    int progen_node;
    bool birth;
};


// To sort the particles by particleID
int struct_cmp(const void *a, const void *b)
{
    struct Frame *ia = (struct Frame *)a;
    struct Frame *ib = (struct Frame *)b;
    return (int) (ia->NodeID - ib->NodeID);
}


// To read halo files
static int display_info(const char *fpath, const struct stat *sb, int tflag, struct FTW *ftwbuf)
{    
    tempLoad = fpath;
    if((tempLoad != "./FileRead.cpp") && (tempLoad != ".") && (tempLoad != "./FileRead") && (tempLoad != "./FileRead.cpp~") && (tempLoad != "./example.txt") && (tempLoad != "./T1T2data_V2.txt") && (tempLoad != "./data_Halo4.dat"))
    {
      listFiles.push_back(tempLoad); cout << "pushed " << endl;   
    }
    cout << "listFiles.size() " << listFiles.size() << endl;
    return 0;           /* To tell nftw() to continue */
}



int main(int argc, char *argv[])
{
    int flags = 0;
    int file_iter, fileDiff = 0, switchStep = 0, transferVect1 = 0, transferVect2 = 0, curFile = 0, ref1 = 0, ref2 = 0, global_comparison = 0, loop_num = 0, wait = 0;
    int Global_halo = 0;
    char * cstr, *token;
    int timeStep;	
    vector<int> VSobj;
    vector<int> VObj;
    vector<int> VNode;
    vector<Halo_param> H_param;
    vector<int> listFilesTime;
    vector<int> lastStep_halos;
    vector<int> death_record;
    vector<int> merger_record;
    vector<int> split_record;
    vector <Merger_info> tempVec;
    vector <int> tempVec_vol;
    vector <float> tempVec_gain;
    vector< vector <Merger_info> > MergerTree;
    vector< vector <int> > MergerTree_vol;
    vector< vector <float> > MergerTree_gain;
    vector< vector <Halo_param> > MergerTree_halo_param;
    vector<int> descendent;
    vector<int> sorted_halo;
    vector<float> galact_time;
    int gt = 0;
    int merger_iter = 0;
    Frame* Ftime1;
    Frame* Ftime2;

    
   if (argc > 2 && strchr(argv[2], 'd') != NULL)
        flags |= FTW_DEPTH;
   
   if (argc > 2 && strchr(argv[2], 'p') != NULL)
        flags |= FTW_PHYS;
    
   if (nftw((argc < 2) ? "." : argv[1], display_info, 20, flags) == -1) 
   {
        perror("nftw");
        exit(EXIT_FAILURE);
    }
    
    cout << "listFiles.size() before sort" << listFiles.size() << endl;
   sort(listFiles.begin(), listFiles.end());
    
    for(int i = 0; i < listFiles.size(); i++)
    {
      cout << "Path in string " << listFiles.at(i)  << " and size " << sizeof(listFiles.at(i)) << endl;
    }
    //exit(EXIT_SUCCESS);
   cout << "listFiles.size() " << listFiles.size() << endl;
   
   
   // To store file names in a vector
    for(int i=0; i< listFiles.size(); i++)
    {
        cstr = new char [listFiles[i].size()+1];
        strcpy (cstr, listFiles[i].c_str());
        token=strtok (cstr,"_");
        token=strtok(NULL,"_");
        timeStep = atoi (token);
        listFilesTime.push_back(timeStep);  
		
		// To store timesteps for writing a tree to file later
		if(i == 0)
		{
			galact_time.push_back(timeStep);
			gt++;
		}
		if(galact_time[gt - 1] != timeStep)
		{
			galact_time.push_back(timeStep);
			gt++;
		}
    }

	
    //Open file to write tracking results for reference
    ofstream myfile;
    myfile.open("T1T2data_V2.txt");


    float pos1=0, vel1=0, pos2=0, vel2=0, pos3=0, vel3=0, tag1=0; //Positions, velocities and subhalo_tags
    int tag2=0; //Particle tag
    float a_pos1 = 0, a_pos2 = 0, a_pos3 = 0, a_vel1 = 0, a_vel2 = 0, a_vel3 = 0; //Average position and velocity of halo
    int increment = 0;
    int high_sub1 = 0, high_sub2 = 0;
    int Halo_count1 = 0, Halo_count2 = 0;
    int countN=0, countO=0, Cnode, Cobj, count1=0, countN2=0, countO2=0, Cnode1, Cobj1;
    int deaths = 0;
    vector<int> HaloVol1;
    vector<int> HaloVol2;
	

    do {      
      
        loop_num++;

        // For the first file of a time step
        if((curFile == 0) || (ref1 == 0))
        {

            ifstream file (listFiles[curFile].c_str(), ios::in|ios::binary);            
            countN = 0;
            int iter = 0; //To categorise the particle parameters according to their position in the file

            if(file.is_open())
            {
     
                while(!file.eof())
                {

                    if ((iter % 8) == 0)
                    {
                        file.read((char *)(&pos1), sizeof(pos1));
                        a_pos1 += pos1;
                    }
     

                    else if ((iter % 8) == 1)
                    {
                        file.read((char *)(&vel1), sizeof(vel1));
                        a_vel1 += vel1;
                    }
 

                    else if ((iter % 8) == 2)
                    {
                        file.read((char *)(&pos2), sizeof(pos2));
                        a_pos2 += pos2;
                    }


                    else if ((iter % 8) == 3)
                    {
                        file.read((char *)(&vel2), sizeof(vel2));
                        a_vel2 += vel2;
                    }


                    else if ((iter % 8) == 4)
                    {
                        file.read((char *)(&pos3), sizeof(pos3));
                        a_pos3 += pos3;
                    }


                    else if ((iter % 8) == 5)
                    {
                        file.read((char *)(&vel3), sizeof(vel3));
                        a_vel3 += vel3;
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
                        if(curFile > 1)
                        {
                            // Use updated halo info from last steps OverlapTable
                            VObj.push_back(lastStep_halos[Halo_count1]); //Assign halo ID
                        }
                        else
                        {
                            VObj.push_back(curFile+1); //Assign halo ID
                        }
                                                //Ftime1[countO].ObjID = filecounter;
                                                //cout << "Subhalo IDs " << VSobj[countO] << "Halo ID " << VObj[countO] << endl;
                        countO++;
                    }
     

                    else if ((iter % 8) == 7) //Particle_tag
                    {
                        file.read((char *)(&tag2), sizeof(tag2));
                        VNode.push_back(tag2);
                                                //Ftime1[countN].NodeID = tag2;
                                                //cout << "Particle ID " << VNode[countN] << endl;
                        countN++;
                    }
                               
                    iter++;

                }

                HaloVol1.push_back(countN);
                Halo_count1++; // Count the number of Halos in time i.
                curFile++; // Keep a count on number of files read.
                ref1++; // The code above was only for the first file.
            }
            file.close();
            file.clear();
           
           
            //calculate the average of all particles for a halo 
            a_pos1 = a_pos1/countN;
            a_pos2 = a_pos2/countN;
            a_pos3 = a_pos3/countN;
            a_vel1 = a_vel1/countN;
            a_vel2 = a_vel2/countN;
            a_vel3 = a_vel3/countN;
                       
         
            increment = (high_sub1 + 1); //Assign value of higest subhalo_tag to continue numbering

            if(curFile == 1)
            {
                Global_halo = Halo_count1;
                tempVec.resize(Global_halo);
                tempVec[Global_halo - 1].host_node = Global_halo;
				tempVec[Global_halo - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
				tempVec[Global_halo - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
            }
           
            H_param.resize(Global_halo);
           
            int push = 0;
            if(curFile>1)
                push = lastStep_halos[Halo_count1 - 1] - 1; //Subtract 1 since Halo_count1 has been incremented and subtract 1 again since we want to put in the elemental position [i-1]
            else
                push = curFile - 1;

			// Store the calculate average position and velocities of the halos
            H_param[push].px = a_pos1;
            H_param[push].py = a_pos2;
            H_param[push].pz = a_pos3;
            H_param[push].vx = a_vel1;
            H_param[push].vy = a_vel2;
            H_param[push].vz = a_vel3;

        }

        else if((switchStep == 0)&&(ref1 == 1))
        {          
            if(listFilesTime[curFile] == listFilesTime[curFile-1])
            {
                ifstream file (listFiles[curFile].c_str(), ios::in|ios::binary);                

                int iter = 0; 
                a_pos1 = 0, a_pos2 = 0, a_pos3 = 0, a_vel1 = 0, a_vel2 = 0, a_vel3 = 0;
                countN = 0;

                if(file.is_open())
                {
     
                    while(!file.eof())
                    {

                        if ((iter % 8) == 0)
                        {
                            file.read((char *)(&pos1), sizeof(pos1));
                            a_pos1 += pos1;
                        }
     

                        else if ((iter % 8) == 1)
                        {
                            file.read((char *)(&vel1), sizeof(vel1));
                            a_vel1 += vel1;
                        }
 

                        else if ((iter % 8) == 2)
                        {
                            file.read((char *)(&pos2), sizeof(pos2));
                            a_pos2 += pos2;
                        }


                        else if ((iter % 8) == 3)
                        {
                            file.read((char *)(&vel2), sizeof(vel2));
                            a_vel2 += vel2;
                        }


                        else if ((iter % 8) == 4)
                        {
                            file.read((char *)(&pos3), sizeof(pos3));
                            a_pos3 += pos3;
                        }


                        else if ((iter % 8) == 5)
                        {
                            file.read((char *)(&vel3), sizeof(vel3));
                            a_vel3 += vel3;
                        }


                        else if ((iter % 8) == 6) //Subhalo_tag
                        {
                            file.read((char *)(&tag1), sizeof(tag1));
                                                                    
                            tag1 = tag1 + increment; // Assigns a global subhalo_tag
                            if(tag1 > high_sub1) // Records the highest subhalo_tag in this file
                            {
                                high_sub1 = tag1;
                            }
                            VSobj.push_back(tag1);
                                                    

                            if(global_comparison > 0)
                            {
                                // Use updated halo info from last step's OverlapTable
                                VObj.push_back(lastStep_halos[Halo_count1]); //Assign halo ID
                            }
                            else if(global_comparison == 0)
                            {
                                VObj.push_back(Global_halo+1); //Assign halo ID                                                    
                            }                          
                            countO++;                          
                        }
     

                        else if ((iter % 8) == 7) //Particle_tag
                        {
                            file.read((char *)(&tag2), sizeof(tag2));
                            VNode.push_back(tag2);
                            countN++;
                        }
                               
                        iter++;

                    }

                    HaloVol1.push_back(countN);
                    Halo_count1++; // Count the number of Halos in time i
                    curFile++; // Keep a count on number of files read
                    if(global_comparison == 0)
                    {
                        Global_halo++; // Keep a track of global halos
                        tempVec.resize(Global_halo);
                        tempVec[Global_halo - 1].host_node = Global_halo;
						tempVec[Global_halo - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
						tempVec[Global_halo - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
                    }
                       
                }
                file.close();
                file.clear();

                //calculate the average of all particles for a halo
                a_pos1 = a_pos1/countN;
                a_pos2 = a_pos2/countN;
                a_pos3 = a_pos3/countN;
                a_vel1 = a_vel1/countN;
                a_vel2 = a_vel2/countN;
                a_vel3 = a_vel3/countN;
               
                increment = (high_sub1 + 1); //Assign value of higest subhalo_tag to continue numbering
                H_param.resize(Global_halo);

                int push = 0;
                if(global_comparison > 0)
                    push = lastStep_halos[Halo_count1 - 1] - 1; //Subtract 1 since Halo_count1 has been incremented and subtract 1 again since we want to put in the elemental position [i-1]
                else
                    push = Global_halo - 1;
               
                H_param[push].px = a_pos1;
                H_param[push].py = a_pos2;
                H_param[push].pz = a_pos3;
                H_param[push].vx = a_vel1;
                H_param[push].vy = a_vel2;
                H_param[push].vz = a_vel3;
               
				if(curFile == listFiles.size())
				{
					MergerTree_halo_param.push_back(H_param);
					switchStep = 3; // Condition met to terminate the do-while loop
				}
					
            }

            else

            {
                switchStep++; // Signal end of files from same timestep
                transferVect1++;
                cout << "TransferVect increased 1" << endl;
            }

        }


        // Transfer particles to the struct      
        if((transferVect1 == 1))
        {
            Ftime1 = new Frame[countO];
          
            for(int sub = 0; sub < countO; sub++) //Trasfer info from vectors
            {
                Ftime1[sub].SObjID = VSobj[sub];
                Ftime1[sub].ObjID = VObj[sub];
                Ftime1[sub].NodeID = VNode[sub];                                                    
            }
            VSobj.clear();
            VObj.clear();
            VNode.clear();

            myfile<<"Ftime1 Counts"<<countO<<" "<<countN<<endl;
 
            //size_t Ftime1_len = sizeof(Ftime1) / sizeof(struct Frame);        

            qsort(Ftime1,countO, sizeof(struct Frame), struct_cmp); // Sort particle list according to particleID
            cout << "Sorted!" << endl;

            if(global_comparison == 0)
            {
                tempVec_vol.resize(Global_halo);
                for(int l = 0; l < HaloVol1.size(); l++)
                    tempVec_vol.at(l) = HaloVol1.at(l);

                MergerTree.push_back(tempVec);
                MergerTree_vol.push_back(tempVec_vol);
                tempVec.clear();
                tempVec_vol.clear();
            }                       

            transferVect1 = 0;

        }



        // No more files from the same time step. Move on to the next time step.
        if((switchStep == 1) && (ref2 == 0))
        {
            increment = 0;
            ifstream file2 (listFiles[curFile].c_str(), ios::in|ios::binary);
            int iter2 = 0; //To categorise the parameters according to their position in the file
            countN2 = 0;

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
                        tag1 = tag1 + increment; // Assigns a global subhalo_tag
                        if(tag1 > high_sub2) // Records the highest subhalo_tag in this file
                        {
                            high_sub2 = tag1;
                        }
                        VSobj.push_back(tag1);                                                                
                        VObj.push_back(Halo_count2+1);
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

                HaloVol2.push_back(countN2);
                Halo_count2++; // Count the number of halos in time i+1.
                curFile++; // Keep a count on number of files read.
                ref2++; // The code above was only for the first file.
            }
            file2.close();
            file2.clear();
     
 
            increment = (high_sub2 + 1); //Assign value of highest subhao_tag to continue numbering
        }

        else if((switchStep == 1) && (ref2 == 1))
        {          
            if((curFile < listFiles.size()) && (listFilesTime[curFile] == listFilesTime[curFile-1]))
            {
                increment = 0;
                ifstream file2 (listFiles[curFile].c_str(), ios::in|ios::binary);
                int iter2 = 0; //To categorise the parameters according to their position in the file
                countN2 = 0;

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
                            tag1 = tag1 + increment; // Assigns a global subhalo_tag
                            if(tag1 > high_sub2) // Records the highest subhalo_tag in this file
                            {
                                high_sub2 = tag1;
                            }
                            VSobj.push_back(tag1);
                            VObj.push_back(Halo_count2+1);
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

                    HaloVol2.push_back(countN2);
                    Halo_count2++; // Count the number of halos in time i+1.
                    curFile++; // Keep a count on number of files read.
                }
                file2.close();
                file2.clear();
     
 
                increment = (high_sub2 + 1); //Assign value of highest subhao_tag to continue numbering
                fileDiff++; //Keep a count of files read in time step i+1

            }

            else
            {
                switchStep++;
                transferVect2++;
                cout << "TransferVect increased 2" << endl;
            }

        }

      
        // Transfer particles to the struct
      
        if(transferVect2 == 1)
        {
            Ftime2 = new Frame[countO2];
            for(int sub2 = 0; sub2 < countO2; sub2++)
            {
                Ftime2[sub2].SObjID = VSobj[sub2];
                Ftime2[sub2].ObjID = VObj[sub2];
                Ftime2[sub2].NodeID = VNode[sub2];
            }
            VSobj.clear();
            VObj.clear();
            VNode.clear();


            myfile<<"Ftime2 Counts"<<countO2<<" "<<countN2<<endl;
            //cout<<"Sizeof"<<(sizeof(Ftime2)/sizeof(Ftime2[0]))<<endl;
            //cout << "Highsub " << high_sub2 << endl; // Comment for faster execution
                      
            size_t Ftime2_len = sizeof(Ftime2) / sizeof(struct Frame);
          

            /*for(int Ncount2 = 0; Ncount2 < countO2; Ncount2++)
            {
                //Halo2.push_back(Ftime2[Ncount2].SObjID);
                HaloVol2[Ftime2[Ncount2].ObjID - 1]++;
            }*/
  

            qsort(Ftime2,countO2, sizeof(struct Frame), struct_cmp);
          
            transferVect2 = 0;
        }
   
      
		// **** Start Overlap Detection ****
        if((switchStep == 2) || (curFile > listFiles.size()))
        {
 
            cout << "deaths " << deaths << " Halo_count1 " << Halo_count1 << " Halo_Count2 " << Halo_count2 << endl;
            vector< vector <int> > OverlapTable ((Halo_count1 + 1 + deaths), vector<int> (Halo_count2 + 1)); // Initialise the Halo Overlap Table

            // Assign dead halo rows their Halo ID's for record continuation
            if (deaths > 0)
            {
                for(int k = 0; k < death_record.size(); k++)
                {                  
                    OverlapTable[death_record[k] - 1][Halo_count2] = death_record[k];                 
                }
            }
            death_record.clear();

            register int iT (0), jT (0);

            // 'While' loop does Overlap detection 
            while (iT<countO || jT<countO2) // Changed from (iT<countO && jT<countO2) so that all particles of both lists are scanned
            {
              
                if ((iT<countO && jT<countO2) && (Ftime1[iT].NodeID == Ftime2[jT].NodeID))  
                {
                  
                    OverlapTable[Ftime1[iT].ObjID - 1][Ftime2[jT].ObjID - 1]++; 
                    OverlapTable[Halo_count1 + deaths][Ftime2[jT].ObjID - 1] = Ftime2[jT].ObjID; // Connect subhalo to its Halo in time step i+1
                    OverlapTable[Ftime1[iT].ObjID - 1][Halo_count2] = Ftime1[iT].ObjID; // Connect subhalo to its Halo in time step i
                  
                    iT++;
                    jT++;
                }

                else
       
                {

                    if((iT >= countO) || (Ftime1[iT].NodeID > Ftime2[jT].NodeID))
                    {
                        OverlapTable[Halo_count1 + deaths][Ftime2[jT].ObjID - 1] = Ftime2[jT].ObjID; // Connect subhalo to its Halo in time step i+1, for subhalo tracking
                                            
                        if (jT < countO2)
                        {                                          
                            jT++;
                        }

                    }
                    else if(Ftime1[iT].NodeID < Ftime2[jT].NodeID)
                    {  
                        OverlapTable[Ftime1[iT].ObjID - 1][Halo_count2] = Ftime1[iT].ObjID; // Connect subhalo to its Halo in time step i, for subhalo tracking

                        if (iT < countO)
                            {
                                iT++;                          
                            }

                    }
                }

            }

			// Thresholding for the special case halos born during split
			for(int iter_r = 0; iter_r < (Halo_count1 + 1 + deaths); iter_r++)
			{
				for(int iter_c = 0; iter_c < (Halo_count2 + 1); iter_c++)
				{
					if((OverlapTable[iter_r][iter_c] > 0) && (OverlapTable[iter_r][iter_c] < 100) && (iter_c != Halo_count2) && (iter_r != (Halo_count1 + deaths)))
					{
						OverlapTable[iter_r][iter_c] = 0;
					}					
				}
			}
			cout << "End of timestep" << endl;
			//cin >> wait;
			
			 
             

          
            // Write tracking results to a file
       
            myfile<< "===========================================================================================" << endl << endl;       
            int new_rows = (Halo_count1 + deaths);
            if(global_comparison == 0)
            {
                tempVec.resize(Global_halo);
                tempVec_vol.resize(Global_halo);
                tempVec_gain.resize(Global_halo);
                merger_record.resize(Global_halo);
				split_record.resize(Global_halo);
            }
          
            deaths = 0;
            float death_merger_add = 0.0;

			// Start scanning for overlap and categorizing events ONE ROW AT A TIME 
            for(int Scount1 = 0; Scount1 < new_rows; Scount1++)
            {              
                vector <int> tempPer;
                vector <int> temp;
                vector <int> tempB;
                vector <int> tempM;
                int checkE = 0, iterB = 0, scanB = 0;
               
                // To check number of overlapping features in timestep i+1
                for(int Scount2 = 0; Scount2 < Halo_count2; Scount2++)
                {
                    if(OverlapTable[Scount1][Scount2]>0)
                    {
						tempPer.push_back(OverlapTable[Scount1][Scount2]);
						temp.push_back(Scount2);
						checkE++; // To check number of overlaps in timestep i+1
                    }
                    
					// Push columns with empty 'first element' into tempB
                    else if((Scount1 == 0) && (OverlapTable[Scount1][Scount2] == 0))
                    {
                        tempB.push_back(Scount2);
                        iterB++;
                    }
                              
                }
                

                   
                float PerOver = 0.0;
                double var1= 0.0, var2 = 0.0;

                // Check for BIRTH event
				// Since this is the first row, if any element's column is fully zero, it is a birth
				if((Scount1 == 0) && (tempB.size() > 0))
                {       
                    for(int birth = 0; birth < tempB.size(); birth++)
                    {
                        scanB = 0; // Reinitialise birth check counter
                        for(int chkB = 0; chkB < new_rows; chkB++)
                        {
                            if(OverlapTable[chkB][tempB[birth]] > 0)
                                scanB++;
                        }
                   
                        if(scanB == 0)
                        {
                            OverlapTable[new_rows][tempB[birth]] = Global_halo + 1; //Use the Global_halo to record birth of a halo in timestep ti+1
                            Global_halo++; // Increment global halo assignment variable
                            tempVec.resize(Global_halo);
                            tempVec_vol.resize(Global_halo);
                            tempVec_gain.resize(Global_halo);
                            merger_record.resize(Global_halo);
							split_record.resize(Global_halo);
                            H_param.resize(Global_halo);
                            myfile << "                                SubHalo ID " << tempB[birth] << " of Halo ID " << OverlapTable[new_rows][tempB[birth]] << "(" << HaloVol2.at(tempB[birth]) << ")" << " is born." << endl << endl;                              
                            tempVec[OverlapTable[new_rows][tempB[birth]] - 1].host_node = OverlapTable[new_rows][tempB[birth]];
							tempVec[OverlapTable[new_rows][tempB[birth]] - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
                            tempVec[OverlapTable[new_rows][tempB[birth]] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
							tempVec[OverlapTable[new_rows][tempB[birth]] - 1].birth = true;
							tempVec_vol.at(OverlapTable[new_rows][tempB[birth]] - 1) = HaloVol2.at(tempB[birth]);
                            tempVec_gain.at(OverlapTable[new_rows][tempB[birth]] - 1) = 0;
                        }
                    }
					tempB.clear(); // Free birth candidates
                }
           


                //DEATH event
                if(checkE == 0)
                {
                    myfile << "SubHalo ID " << Scount1 << " of Halo ID " << OverlapTable[Scount1][Halo_count2] << " dies" << endl;
                    deaths++;
                    tempVec[OverlapTable[Scount1][Halo_count2] - 1].host_node = 0;
					tempVec[OverlapTable[Scount1][Halo_count2] - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
                    tempVec[OverlapTable[Scount1][Halo_count2] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
					tempVec[OverlapTable[Scount1][Halo_count2] - 1].birth = false;
					tempVec_vol.at(OverlapTable[Scount1][Halo_count2] - 1) = 0;
                    tempVec_gain.at(OverlapTable[Scount1][Halo_count2] - 1) = 0;
                    death_record.push_back(OverlapTable[Scount1][Halo_count2]);
                      
                }



                //SPLIT event
                else if((checkE > 1) && ((split_record.at(OverlapTable[Scount1][Halo_count2] - 1) != 1)))
                {
                    myfile<<"SubHalo ID  "<<Scount1<<" ("<</*HaloVol1[Scount1]<<*/") of Halo ID " << OverlapTable[Scount1][Halo_count2] << endl; //NormObj: (Scount1) else use Halo1[Scount1]
                                              
                    vector<int> frag(temp.size());                      
                    frag.at(0)=temp.at(0);    

					// Find biggest fragment
                    for(int split = 0; split < checkE - 1; split++)
                    {
                        if(OverlapTable[Scount1][temp[split + 1]] > OverlapTable[Scount1][frag[0]])
                        {  
                            frag.at(split + 1) = frag.at(0);                               
                            frag.at(0) = temp.at(split + 1);                               
                        }
                        else
                            frag.at(split + 1) = temp.at(split + 1);
                    }
                      

					// Assign haloID's to fragments and update tree
                    for(int split = 0; split < checkE; split++)
                    {
                        var1 = OverlapTable[Scount1][frag[split]];
                        var2 = HaloVol2.at(frag[split]);
                        PerOver = (var1/var2)*100;
                        if(split == 0)
                        {
                            OverlapTable[new_rows][frag[0]] = OverlapTable[Scount1][Halo_count2]; //Biggest fragment continues as the original halo
							tempVec[OverlapTable[new_rows][frag[split]] - 1].birth = false;
                        }
                        else
                        {
                            OverlapTable[new_rows][frag[split]] = Global_halo + 1; //Use the Global_halo to record birth of a halo due to split in timestep ti+1
                            Global_halo++; //Increment
                            tempVec.resize(Global_halo);
                            tempVec_vol.resize(Global_halo);
                            tempVec_gain.resize(Global_halo);
                            merger_record.resize(Global_halo);
							split_record.resize(Global_halo);
                            H_param.resize(Global_halo);
							tempVec[OverlapTable[new_rows][frag[split]] - 1].birth = true;
                        }

                        tempVec[OverlapTable[new_rows][frag[split]] - 1].host_node = OverlapTable[new_rows][frag[split]];
						tempVec[OverlapTable[new_rows][frag[split]] - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
                        tempVec[OverlapTable[new_rows][frag[split]] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module							
						tempVec_vol.at(OverlapTable[new_rows][frag[split]] - 1) = HaloVol2.at(frag[split]);
                        tempVec_gain.at(OverlapTable[new_rows][frag[split]] - 1) = PerOver;
						split_record.at(OverlapTable[Scount1][Halo_count2] - 1) = 1;
                        myfile<<"                     ->> ("<<tempVec_gain.at(OverlapTable[new_rows][frag[split]] - 1)<<"%) ->> SubHalo ID "<<frag[split]<<" ("<< HaloVol2.at(frag[split]) <<") of Halo ID " << OverlapTable[new_rows][frag[split]] <<endl; //NormObj: temp[split] else use Halo2[Scount2]
                        //cout<<"Split"<<endl;
                          
                    }


					// Check for Split-Merge event
                    for(int split = 0; split < temp.size(); split++)
                    {
						tempM.clear();
						tempM.push_back(OverlapTable[new_rows][temp[split]]);
                        int countM = 0;		

						// Scan rows to check if merging elements present in column
                        for(int merge = 0; merge < new_rows; merge++)
                        {
                            if(OverlapTable[merge][temp[split]] > 0)
                            {
                                if(merge != Scount1)
									tempM.push_back(OverlapTable[merge][Halo_count2]);
								countM++;                                                                   
                            }
                        }


                        //Merge Part
                        if((countM > 1) && (merger_record.at(OverlapTable[Scount1][Halo_count2] - 1) != 1))
                        {
							for(int show = 0; show < tempM.size(); show++)
							{
								cout << " tempM = " << tempM.at(show) << endl;
							}

                            int TotOvHalo = 0, TotHalo1 = 0;
                            vector<int> merge_frag(tempM.size()); //Just check this tempM.size if debugging fails
                            merge_frag.at(0) = tempM.at(0);
							
							// Store smaller halo ID in merge_frag[0]
							for(int Cmerge = 0; Cmerge < tempM.size() - 1; Cmerge++)
							{                               
								if(tempM[Cmerge + 1] < merge_frag[0])
								{
									merge_frag.at(Cmerge + 1) = merge_frag.at(0);                                
									merge_frag.at(0) = tempM.at(Cmerge + 1);
								}
								else
								{
									merge_frag.at(Cmerge + 1) = tempM.at(Cmerge + 1);
								}                               
							}

								
							// Assign haloID's to merging halos
							for(int Cmerge = 0; Cmerge < tempM.size(); Cmerge++)
                            {
                                float PerInOver = 0.0;
                                double InVar1 = 0.0, InVar2 = 0.0;
                                                                   
                                if(Cmerge == 0)
                                {
									TotOvHalo = OverlapTable[Scount1][temp[split]] + TotOvHalo;
									PerInOver = (InVar1/InVar2)*100;
										
									tempVec[tempM[Cmerge] - 1].host_node = merge_frag[0];
									tempVec[tempM[Cmerge] - 1].merge_node = merge_frag[Cmerge + 1];
									tempVec[tempM[Cmerge] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
									tempVec[tempM[Cmerge] - 1].birth = false;
                                    tempVec_vol.at(merge_frag[Cmerge] - 1) = HaloVol2.at(temp[split]); // Min halo ID retains all particles
                                    tempVec_gain.at(merge_frag[Cmerge] - 1) = PerInOver;
                                    myfile<<"SubHalo ID "<<tempM[Cmerge]<<" ("<</*HaloVol1[tempM[Cmerge]]<<*/ ") of Halo ID " << merge_frag[0] << " (" << /*tempVec_gain.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) <<*/ "%) " <<endl; //NormObj: tempM[Cmerge] else use Halo1[tempM[Cmerge]]                                   																			
								}
                                else
                                {									
									TotOvHalo = OverlapTable[tempM[Cmerge]][temp[split]] + TotOvHalo;
									PerInOver = (InVar1/InVar2)*100;										
									tempVec[tempM[Cmerge] - 1].host_node = merge_frag[0];
									tempVec[tempM[Cmerge] - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
									tempVec[tempM[Cmerge] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
									tempVec[tempM[Cmerge] - 1].birth = false;
                                    tempVec_vol.at(merge_frag[Cmerge] - 1) = 0; // Other participants die										
                                    tempVec_gain.at(merge_frag[Cmerge] - 1) = PerInOver;	
									merger_record.at(tempM[Cmerge] - 1) = 1; // Record merger for participants
									cout << "inside else, temp M = " << tempM[Cmerge] << " temp[split] = " << temp[split] << endl;
                                    myfile<<"SubHalo ID "<<tempM[Cmerge]<<" ("<</*HaloVol1[tempM[Cmerge]]<<*/ ") of Halo ID " << merge_frag[Cmerge] << " (" << /*tempVec_gain.at(OverlapTable[new_rows][temp[split]] - 1) <<*/ "%) " <<endl; //NormObj: tempM[Cmerge] else use Halo1[tempM[Cmerge]]                                   
                                    death_record.push_back(merge_frag[Cmerge]);        										
									deaths++;
                                }//cout << "In var1 " << Cmerge << " tempVec " << tempM.size() << endl;cin >> wait;                  
																		
                            }						
                                
                            var1 = TotOvHalo;
                            var2 = HaloVol2.at(temp[split]);	
                            PerOver = (var1/var2)*100;
                            OverlapTable[new_rows][temp[split]] = merge_frag[0]; // Give the min halo ID to the column (next timeste halo) 
							myfile<<"                                 }-> ("<<PerOver<<"%) }-> SubHalo ID "<<temp[split]<<" ("<</*HaloVol2[temp[0]]<<*/") of Halo ID " << OverlapTable[new_rows][temp[split]] <<endl; //NormObj: temp[0] else use Halo2[temp[0]]			
                           
                            merge_frag.clear();
							merger_iter++;
                        }
							
                    }

                    frag.clear();
                }



                // *** Condition for MERGE/MERGE-SPLIT/CONTINUATION event ***
                else if((checkE == 1) && (merger_record.at(OverlapTable[Scount1][Halo_count2] - 1) != 1))
                {
                    int countM = 0;
					vector<int> tempS;
					vector<int> rec_merge_split;

					// Store rows that have non-zero element
                    for(int merge = 0; merge < new_rows; merge++)
                    {
                        if(OverlapTable[merge][temp[0]]>0)
                        {
                            tempM.push_back(merge);
                            countM++;                              
                        }
                    }


					// *** MERGE-SPLIT EVENT ***
					for(int i = 0; i < tempM.size(); i++)
					{
						int countS = 0;

						// Find if the row has splitting halos by scanning each element of the row
						for(int j = 0; j < Halo_count2; j++)
						{
							if(OverlapTable[tempM[i]][j]>0)
                            {
                                tempS.push_back(j);
								countS++;
                                                                   
                            }								
						}

						// If splitting halos found, sort halos according to size and name them
						if(countS > 1) 
						{
							rec_merge_split.push_back(tempM[i]); // Record the rows that have splitting halos
							myfile<<"SubHalo ID  "<<tempM[i]<<" ("<</*HaloVol1[Scount1]<<*/") of Halo ID " << OverlapTable[tempM[i]][Halo_count2] << endl; //NormObj: (Scount1) else use Halo1[Scount1]
							vector<int> frag(tempS.size());						
							frag.at(0)=tempS.at(0);

							for(int split = 0; split < tempS.size() - 1; split++)
							{
								if(OverlapTable[tempM[i]][tempS[split + 1]] > OverlapTable[tempM[i]][frag[0]])
								{	
									frag.at(split + 1) = frag.at(0); 								
									frag.at(0) = tempS.at(split + 1); 								
								}
								else
									frag.at(split + 1) = tempS.at(split + 1);
							}


							for(int split = 0; split < tempS.size(); split++)
							{
								var1 = OverlapTable[tempM[i]][frag[split]];									
								var2 = HaloVol2.at(frag[split]);
								PerOver = (var1/var2)*100;
								if(split == 0)
								{
									OverlapTable[new_rows][frag[0]] = OverlapTable[tempM[i]][Halo_count2]; //Biggest fragment continues as the original halo
									tempVec[OverlapTable[new_rows][frag[split]] - 1].birth = false;
								}
								else
								{
									OverlapTable[new_rows][frag[split]] = Global_halo + 1; //Use the Global_halo to record birth of a halo due to split in timestep ti+1
									Global_halo++; //Increment
									tempVec.resize(Global_halo);
									tempVec_vol.resize(Global_halo);
									tempVec_gain.resize(Global_halo);
									merger_record.resize(Global_halo);
									split_record.resize(Global_halo);
									H_param.resize(Global_halo);
									tempVec[OverlapTable[new_rows][frag[split]] - 1].birth = true;
								}
																		
								tempVec[OverlapTable[new_rows][frag[split]] - 1].host_node = OverlapTable[new_rows][frag[split]];
								tempVec[OverlapTable[new_rows][frag[split]] - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
								tempVec[OverlapTable[new_rows][frag[split]] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
								tempVec_vol.at(OverlapTable[new_rows][frag[split]] - 1) = HaloVol2.at(frag[split]);
								tempVec_gain.at(OverlapTable[new_rows][frag[split]] - 1) = PerOver;
								split_record.at(OverlapTable[tempM[i]][Halo_count2] - 1) = 1;
								myfile<<"                     ->> ("<<tempVec_gain.at(OverlapTable[new_rows][frag[split]] - 1)<<"%) ->> SubHalo ID "<<frag[split]<<" ("<< HaloVol2.at(frag[split]) <<") of Halo ID " << OverlapTable[new_rows][frag[split]] <<endl; //NormObj: temp[split] else use Halo2[Scount2]							
								//cout<<"Split"<<endl;
							
							}
							frag.clear();
						}	
						tempS.clear();
					}
						


					// After 'splitting halos' are named, enter the merger module taking into consideration splitting and merging halos
					if((countM > 1) && (rec_merge_split.size() > 0))
					{							
						tempS.push_back(OverlapTable[new_rows][temp[0]]); // Push the column where we detected a single filled element
						int countM = 0;							
						for(int merge = 0; merge < tempM.size(); merge++)
						{
							int split_check = 0; 
							for(int d = 0; d < rec_merge_split.size(); d++) // Do not include the row which has splitting halos in the Min HaloID calculations
							{
								if(tempM[merge] == rec_merge_split[d])
								{
									split_check++;
								}
							}									

								if(split_check == 0)
									tempS.push_back(OverlapTable[tempM[merge]][Halo_count2]);
								countM++;                                         									
						}	

							
						vector<int> merge_frag(tempS.size());
						merge_frag.at(0) = tempS.at(0);

						// Find Min HaloID
						for(int Cmerge = 0; Cmerge < tempS.size() - 1; Cmerge++)
						{
							if(tempS[Cmerge + 1] < merge_frag[0])
							{
								merge_frag.at(Cmerge + 1) = merge_frag.at(0);
								merge_frag.at(0) = tempS.at(Cmerge + 1);
							}
							else
							{
								merge_frag.at(Cmerge + 1) = tempS.at(Cmerge + 1);
							}
						}

						// Relabel halos in next timestep
						int TotOvHalo = 0, TotHalo1 = 0;
						for(int Cmerge = 0; Cmerge < tempS.size(); Cmerge++)
						{
							float PerInOver = 0.0;
							double InVar1 = 0.0, InVar2 = 0.0;
							InVar1 = OverlapTable[tempM[Cmerge]][temp[0]];
							InVar2 = HaloVol2.at(temp[0]);
							PerInOver = (InVar1/InVar2)*100;

							if(Cmerge == 0)
							{
								TotOvHalo = OverlapTable[tempM[Cmerge]][temp[0]] + TotOvHalo;
								tempVec[tempS[Cmerge] - 1].host_node = merge_frag[0];
                                tempVec[tempS[Cmerge] - 1].merge_node = merge_frag[Cmerge + 1];
                                tempVec[tempS[Cmerge] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
								tempVec[tempS[Cmerge] - 1].birth = false;
								tempVec_vol.at(merge_frag[0] - 1) = HaloVol2.at(temp[0]);
								tempVec_gain.at(merge_frag[0] - 1) = PerInOver;
								merger_record.at(tempM[Cmerge]) = 1; // Since here tempM stores row number
								myfile<< "SubHalo ID " << tempM[Cmerge] << " (" << /*HaloVol1[tempS[Cmerge]] <<*/ ") of HaloID" << merge_frag[0] << " (" << /*tempVec_gain.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) <<*/ "%) " << endl;
							}
							else
							{
								TotOvHalo = OverlapTable[tempM[Cmerge]][temp[0]] + TotOvHalo;
								tempVec[tempS[Cmerge] - 1].host_node = merge_frag[0];
								tempVec[tempS[Cmerge] - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
                                tempVec[tempS[Cmerge] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
								tempVec[tempS[Cmerge] - 1].birth = false;
								tempVec_vol.at(merge_frag[Cmerge] - 1) = 0;
								tempVec_gain.at(merge_frag[Cmerge] - 1) = 0;
								merger_record.at(tempM[Cmerge]) = 1; // Since here tempM stores row number
								myfile << "SubHalo ID " << tempM[Cmerge] << " (" << /*HaloVol1[tempS[Cmerge]] <<*/ ") of HaloID" << merge_frag[Cmerge] << " (" << /*tempVec_gain.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) <<*/ "%) " << endl;
								death_record.push_back(merge_frag[Cmerge]);
								deaths++;
							}
							merger_iter++;
						}

						var1 = TotOvHalo;
                        var2 = HaloVol2.at(temp[0]);
						PerOver = (var1 / var2) * 100;
						OverlapTable[new_rows][temp[0]] = merge_frag[0];
						myfile<<"                                 }-> (" << PerOver << "%) }-> SubHalo ID " << Scount1 << " (" <</*HaloVol2[temp[0]]<<*/") of Halo ID " << OverlapTable[new_rows][temp[0]] <<endl;
						merge_frag.clear();
					}
					tempS.clear();
					rec_merge_split.clear();


                    // *** ONLY-MERGE EVENT ***
                    if((countM>1) && (merger_record.at(OverlapTable[Scount1][Halo_count2] - 1) != 1))
                    {
                        int TotOvHalo = 0, TotHalo1 = 0;
                        vector<int> merge_frag(tempM.size()); //Just check this tempM.size if debugging fails
                        merge_frag.at(0) = tempM.at(0);
                           

                        for(int Cmerge = 0; Cmerge < countM - 1; Cmerge++)
                        {              
                              
                            if(OverlapTable[tempM[Cmerge + 1]][Halo_count2] < OverlapTable[merge_frag[0]][Halo_count2])
                            {
                                merge_frag.at(Cmerge + 1) = merge_frag.at(0);                               
                                merge_frag.at(0) = tempM.at(Cmerge + 1);
                            }
                            else
                            {
                                merge_frag.at(Cmerge + 1) = tempM.at(Cmerge + 1);
                            }
                              
                        }


                        for(int Cmerge = 0; Cmerge < countM; Cmerge++)
                        {
                            float PerInOver = 0.0;
                            double InVar1 = 0.0, InVar2 = 0.0;
                            InVar1 = OverlapTable[tempM[Cmerge]][temp[0]];
                            InVar2 = HaloVol2.at(temp[0]);
                            PerInOver = (InVar1/InVar2)*100;
                              
                            if(Cmerge == 0)
                            {
                                tempVec[OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1].host_node = OverlapTable[merge_frag[0]][Halo_count2];
                                tempVec[OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1].merge_node = OverlapTable[merge_frag[Cmerge + 1]][Halo_count2];
                                tempVec[OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
								tempVec[OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1].birth = false;
								tempVec_vol.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) = HaloVol2.at(temp[0]);
                                tempVec_gain.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) = PerInOver;
                                myfile<<"SubHalo ID "<<tempM[Cmerge]<<" ("<</*HaloVol1[tempM[Cmerge]]<<*/ ") of Halo ID " << OverlapTable[merge_frag[Cmerge]][Halo_count2] << " (" << tempVec_gain.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) << "%) " <<endl; //NormObj: tempM[Cmerge] else use Halo1[tempM[Cmerge]]                 
                            }
                            else
                            { 
                                tempVec[OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1].host_node = OverlapTable[merge_frag[0]][Halo_count2];
								tempVec[OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
                                tempVec[OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
								tempVec[OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1].birth = false;
								tempVec_vol.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) = 0;
                                tempVec_gain.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) = 0;
                                myfile<<"SubHalo ID "<<tempM[Cmerge]<<" ("<</*HaloVol1[tempM[Cmerge]]<<*/ ") of Halo ID " << OverlapTable[merge_frag[Cmerge]][Halo_count2] << " (" << tempVec_gain.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) << "%) " <<endl; //NormObj: tempM[Cmerge] else use Halo1[tempM[Cmerge]]
                                death_record.push_back(OverlapTable[merge_frag[Cmerge]][Halo_count2]);   
                                merger_record.at(OverlapTable[merge_frag[Cmerge]][Halo_count2] - 1) = 1;
                                deaths++;
                                  
                            }
                              

                            TotOvHalo = OverlapTable[tempM[Cmerge]][temp[0]] + TotOvHalo; //For percentage calculations
                            /*TotHalo1 = HaloVol1[tempM[Cmerge]] + TotHalo1;*/ //For percentage calculations
                        }


                        var1 = TotOvHalo;
                        var2 = TotHalo1;
                        PerOver = (var1/var2)*100;
                        OverlapTable[new_rows][temp[0]] = OverlapTable[merge_frag[0]][Halo_count2];                   
                        myfile<<"                                 }-> ("<<PerOver<<"%) }-> SubHalo ID "<<temp[0]<<" ("<</*HaloVol2[temp[0]]<<*/") of Halo ID " << OverlapTable[new_rows][temp[0]] <<endl; //NormObj: temp[0] else use Halo2[temp[0]]
                           
                        merge_frag.clear();
						merger_iter++;
                    }


                    // *** CONTINUATION EVENT ***
                    else if((countM == 1) && (merger_record.at(OverlapTable[Scount1][Halo_count2] - 1) != 1))
                    {
                        var1 = OverlapTable[Scount1][temp[0]];
                        var2 = HaloVol2.at(temp[0]);
                        PerOver = (var1/var2)*100; //Percentage overlap
                        OverlapTable[new_rows][temp[0]] = OverlapTable[tempM[0]][Halo_count2]; //Assign HaloID of timestep ti to the continuing halo in ti+1
                        tempVec[OverlapTable[new_rows][temp[0]] - 1].host_node = OverlapTable[new_rows][temp[0]];
						tempVec[OverlapTable[new_rows][temp[0]] - 1].merge_node = 0; // Added to avoid error while reading .merge_node value in concatinated tree module
                        tempVec[OverlapTable[new_rows][temp[0]] - 1].progen_node = 0; // Added to assign progenitor number to descendents in GALACTICUS module
						tempVec[OverlapTable[new_rows][temp[0]] - 1].birth = false;
						tempVec_vol.at(OverlapTable[new_rows][temp[0]] - 1) = HaloVol2.at(temp[0]);
                        tempVec_gain.at(OverlapTable[new_rows][temp[0]] - 1) = PerOver;
                        myfile<<"SubHalo ID "<<Scount1<<" ("<</*HaloVol1[Scount1]<<*/") of Halo ID " << OverlapTable[Scount1][Halo_count2] << " -> ("<<tempVec_gain.at(OverlapTable[new_rows][temp[0]] - 1)<<"%) -> SubHalo ID "<<temp[0]<<" ("<<HaloVol2.at(temp[0])<<") of Halo ID " << OverlapTable[new_rows][temp[0]] <<endl; //NormObj: Scount1, temp[0]
                    //Changed Halo1[Scount1] to Scount1
                    }
                }


                myfile << endl<<endl;
                temp.clear();
                tempPer.clear();                
                tempM.clear();
            }


            MergerTree.push_back(tempVec);
            MergerTree_vol.push_back(tempVec_vol);
            MergerTree_gain.push_back(tempVec_gain);
            MergerTree_halo_param.push_back(H_param);


            lastStep_halos.clear();
            // Save info about unique halo naming to be used in next comparison
            for(int i = 0; i < Halo_count2; i++)
            {
                lastStep_halos.push_back(OverlapTable[new_rows][i]);
                cout << "To next step " << OverlapTable[new_rows][i] << endl;
                                  
            }

            HaloVol1.clear();
            HaloVol2.clear();
            OverlapTable.clear();
            merger_record.clear();
			split_record.clear();
            H_param.clear();
            tempVec.clear();

            ref1 = 0; // Reinitialize for new comparison iteration
            ref2 = 0;
            high_sub1 = 0; high_sub2 = 0;
            Halo_count1 = 0; Halo_count2 = 0;
            countN=0, countO=0, Cnode, Cobj, count1=0, countN2=0, countO2=0, Cnode1, Cobj1, wait;

            // Conditions to either continue or terminate code
            switchStep = 0; // Start with next iteration 
            curFile = curFile - fileDiff - 1;
            fileDiff = 0;
            global_comparison++;
        }


    } while(switchStep != 3);


    //Print the Merger Tree
    myfile << endl << endl;
	cout << "Done tracking and number of mergers " << merger_iter << endl;
    /*for(int y = 0; y < MergerTree.size()-1; y++)
    {
        for(int x = 0; x < MergerTree[y].size(); x++)
        {
            myfile << MergerTree[y][x].host_node << endl;
            myfile << "(" << MergerTree_vol[y][x] << ")" << endl;
            //if(y < MergerTree.size() - 1)
                myfile << MergerTree_halo_param[y][x].px << endl << MergerTree_halo_param[y][x].py << endl << MergerTree_halo_param[y][x].pz << endl << MergerTree_halo_param[y][x].vx << endl << MergerTree_halo_param[y][x].vy << endl << MergerTree_halo_param[y][x].vz << endl << endl;
        }
        myfile << endl <<"*****************************"<< endl;
    }*/

    // To pull halo history with max mass
    /*int max_vol = MergerTree_vol[MergerTree_vol.size() - 2][0]; // - 2 because MergerTree_vol and MergerTree start from 88-165 as opposed to MergerTree_halo_param: 88-155 
    int max_halo = 0;
    for(int y = 0; y < MergerTree_vol[MergerTree_vol.size() - 1].size(); y++)
    {
        if(MergerTree_vol[MergerTree_vol.size() - 1][y] > max_vol)
        {
            max_vol = MergerTree_vol[MergerTree_vol.size() - 1][y];
            max_halo = y;           
        }
    }*/
	cout << "Max mass halo is " << /*max_halo << " with mass " << max_vol <<*/ endl;
    cin >> wait;
 


	// Assigning serial halo number and descendent number to halos in merger tree for GALACTICUS
	int g = 0; // Uncomment for building a merger tree
	int sorted_edge = 120; // An increment counter to sort 
	vector< vector <float> > Galact; // Vector housing serial halo number and descendent number
	vector <float> Galact_row; // Temp vector to push every halos information into Galact

	//for(int fst = 2; fst < MergerTree_vol[MergerTree_vol.size() - 2].size(); fst++) // FOR loop is to push trees of all halos present at last step
	//{//cout << " y is " << fst << endl;
		//if(MergerTree_vol[MergerTree_vol.size() - 1][fst] > 0) // IF loop detects halo with mass>0 (if halo exists at last step)
		//{
			for(int x = MergerTree.size() - 1; x > -1; x--)
			{
				if(x == MergerTree.size() - 1) // For the first/ROOT halo
				{
					if(MergerTree[x][fst].host_node != 0)
					{
						descendent.push_back(MergerTree[x][fst].host_node);
						sorted_halo.push_back(sorted_edge);
						sorted_edge = sorted_edge - 10;
						MergerTree[x - 1][MergerTree[x][fst].host_node - 1].progen_node = g;				
					}

					if(MergerTree[x][fst].merge_node > 0)
					{			
						MergerTree[x - 1][MergerTree[x][fst].merge_node - 1].progen_node = g;
					}
					g++;
			
					Galact_row.push_back(MergerTree[x - 1][MergerTree[x][fst].host_node - 1].progen_node); // HaloID
					Galact_row.push_back(MergerTree[x][fst].host_node); // Original HaloID
					Galact_row.push_back(-1); // Progen HaloID
					Galact_row.push_back((1/(0.0196078 + (galact_time[x] * (1.96 / 1000)))) - 1); 
					Galact_row.push_back(MergerTree_vol[x][fst]);
					Galact_row.push_back(MergerTree_halo_param[x][fst].px);
					Galact_row.push_back(MergerTree_halo_param[x][fst].py);
					Galact_row.push_back(MergerTree_halo_param[x][fst].pz);
					Galact_row.push_back(MergerTree_halo_param[x][fst].vx);
					Galact_row.push_back(MergerTree_halo_param[x][fst].vy);
					Galact_row.push_back(MergerTree_halo_param[x][fst].vz);cout << Galact_row[0] << " (" << Galact_row[1] << ") " << "*" << Galact_row[2] << "  T" << Galact_row[3] << " " << Galact_row[4] << " " << Galact_row[5] << " " << Galact_row[6] << " " << Galact_row[7] << " " << Galact_row[8] << " " << Galact_row[9] << " " << Galact_row[10] << " ";
            
					Galact.push_back(Galact_row);
					Galact_row.clear();
				}
				else
				{
					for(int n = 0; n < descendent.size(); n++)
					{
						if(MergerTree[x + 1][descendent.at(n) - 1].merge_node > 0)
						{
							descendent.push_back(MergerTree[x + 1][descendent.at(n) - 1].merge_node);
							sorted_halo.push_back(sorted_edge);
							sorted_edge = sorted_edge - 10;
						}
					}
			 
					for(int n = 0; n < descendent.size(); n++)
					{
						if(MergerTree[x][descendent.at(n) - 1].host_node != 0)
						{					
							Galact_row.push_back(g); // HaloID
							Galact_row.push_back(MergerTree[x][descendent.at(n) - 1].host_node); // Original HaloID
							Galact_row.push_back(MergerTree[x][descendent.at(n) - 1].progen_node); // Progen HaloID													
							Galact_row.push_back((1/(0.0196078 + (galact_time[x] * (1.96 / 1000)))) - 1);
							Galact_row.push_back(MergerTree_vol[x][descendent.at(n) - 1]);
							Galact_row.push_back(MergerTree_halo_param[x][descendent.at(n) - 1].px);
							Galact_row.push_back(MergerTree_halo_param[x][descendent.at(n) - 1].py);
							Galact_row.push_back(MergerTree_halo_param[x][descendent.at(n) - 1].pz);
							Galact_row.push_back(MergerTree_halo_param[x][descendent.at(n) - 1].vx);
							Galact_row.push_back(MergerTree_halo_param[x][descendent.at(n) - 1].vy);
							Galact_row.push_back(MergerTree_halo_param[x][descendent.at(n) - 1].vz);cout << Galact_row[0] << " (" << Galact_row[1] << ") " << "*" << Galact_row[2] << "  T" << Galact_row[3] << " " << Galact_row[4] << " " << Galact_row[5] << " " << Galact_row[6] << " " << Galact_row[7] << " " << Galact_row[8] << " " << Galact_row[9] << " " << Galact_row[10] << " ";
							Galact.push_back(Galact_row);
							g++;
						}
				
						if((MergerTree[x][descendent.at(n) - 1].merge_node > 0) && (x > 0)) // If halo has a merging part in earlier timestep
						{					
							MergerTree[x - 1][descendent.at(n) - 1].progen_node = Galact_row.at(0);
							MergerTree[x - 1][MergerTree[x][descendent.at(n) - 1].merge_node - 1].progen_node = Galact_row.at(0);
						}
						else if (MergerTree[x][descendent.at(n) - 1].birth != true) // If halo exits in earlier timestep
						{
							if((MergerTree[x - 1][descendent.at(n) - 1].host_node > 0) && (x > 0))
								MergerTree[x - 1][descendent.at(n) - 1].progen_node = Galact_row.at(0);
						}
						else if((MergerTree[x][descendent.at(n) - 1].birth == true) && (x > 0)) // If halo was born in this timestep
						{
							descendent.erase(descendent.begin() + n);
							sorted_halo.erase(sorted_halo.begin() + n);
							n--;
						}					
				
						Galact_row.clear();		
				
					}
				} cout << endl << endl;
			}//cout << " Before Writing file " << endl;
			
		//}descendent.clear();
	//}




	// *********** Write to a binary file*******
	ofstream os ("data_Halo4.dat", ios::binary);
	for(int r = 0; r < Galact.size(); r++)
	{
		int Gsize = Galact[r].size(); cout << "Galact[1].size()" << Galact[r].size() << endl;
		os.write((const char*)&Galact[r][0], sizeof(Galact[r][0]) * (Gsize));
	}	
	os.close();


	// **** Read file ****
	/*float bin;
	ifstream readfile ("data_Halo32.dat", ios::in|ios::binary); // Reading the file shows that the last variable is getting stored/printed twice. Look that up.
	if(readfile.is_open())
	{
		while(!readfile.eof())
		{
			readfile.read((char *)(&bin), sizeof(bin));
			cout << "Bin " << bin << endl;
			cin >> wait;
		}
	
	
	}*/
   
   
   
}