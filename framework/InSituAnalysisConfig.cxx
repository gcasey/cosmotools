#include "CosmologyToolsConfiguration.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

namespace cosmotk
{

CosmologyToolsConfiguration::CosmologyToolsConfiguration()
{
  this->Communicator  = MPI_COMM_NULL;
  this->RawConfigData = "";
  this->ConfigFile    = "";
  this->Version       = 0.0;
  this->Visualization = false;
  this->Server        = "127.0.0.1";
  this->Port          = 2222;
  this->VisualizationFrequency = 20;
}

//------------------------------------------------------------------------------
CosmologyToolsConfiguration::~CosmologyToolsConfiguration()
{
  this->Clear();
}

//------------------------------------------------------------------------------
void CosmologyToolsConfiguration::Clear()
{
  this->RawConfigData = "";
  this->ToolNames.clear();
  this->ToolStatus.clear();
  this->ToolToDictionary.clear();
}

//------------------------------------------------------------------------------
void CosmologyToolsConfiguration::ReadInConfigurationFile()
{
  std::stringstream sstream;

  std::ifstream ifs;
  ifs.open(this->ConfigFile.c_str());
  if( !ifs.is_open() )
    {
    std::cerr << "ERROR: Cannot open file: " << this->ConfigFile << std::endl;
    std::cerr << "FILE: " << __FILE__ << std::endl;
    std::cerr << "LINE: " << __LINE__ << std::endl;
    MPI_Abort(this->Communicator,-1);
    } // END if

  // Read all the contents of the file in a stringstream
  sstream.str(std::string(""));
  sstream << ifs.rdbuf();
  ifs.close();

  this->RawConfigData = sstream.str();
}

//------------------------------------------------------------------------------
void CosmologyToolsConfiguration::RegisterAnalysisTool(
        std::string toolName, bool status)
{
  this->ToolNames.push_back( toolName );
  this->ToolStatus[ toolName ] = status;

  Dictionary toolDictionary;
  this->ToolToDictionary[ toolName ] = toolDictionary;
}

//------------------------------------------------------------------------------
bool CosmologyToolsConfiguration::Parsebool(std::string s)
{
  assert("pre: string parameter is empty" && (!s.empty()));
  std::transform(s.begin(),s.end(),s.begin(),::toupper);
  if( (s == "YES") || (s == "TRUE") || (s == "ON") || (s=="ENABLED") )
    {
    return true;
    }
  return false;
}

//------------------------------------------------------------------------------
bool CosmologyToolsConfiguration::ParseRawData()
{
  assert("pre: RawData is empty!" && (!this->RawConfigData.empty()) );

  std::stringstream parser(this->RawConfigData);

  std::string line = "";
  int counter = 0;

  bool detectedVersion = false;
  std::string currentAnalysisTool = "";
  while(std::getline(parser,line))
    {
    ++counter;

    // STEP 0: Tokenize line
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter< std::vector<std::string> >(tokens));
    //BEGIN DEBUG
//        std::cout << counter << ": " << line << std::endl;
//        std::cout << counter << ": NUMTOKENS: " << tokens.size() << std::endl;
//        for(size_t i=0; i < tokens.size(); ++i )
//          {
//          std::cout << "\t" << tokens[i] << std::endl;
//          } // END for all tokens
    //END DEBUG

    // STEP 1: Parse tokens
    if( tokens.size() == 0 )
      {
      continue;
      }
    else if(tokens[0] == "VERSION" )
      {
      assert("pre: cannot parse VERSION tag!" && tokens.size()==2);
      assert("pre: empty VERSION tag!" && !tokens[1].empty());

      detectedVersion = true;
      this->Version = std::atof(tokens[1].c_str());
      }
    else if(tokens[0] == "VISUALIZATION")
      {
      assert("pre: cannot parse VISUALIZATION tag!" && tokens.size()==2);
      assert("pre: empty VISUALIZATION tag!" && !tokens[1].empty());

      this->Visualization = this->Parsebool( tokens[1] );
      }
    else if(tokens[0] == "VIZ_SERVER")
      {
      assert("pre: cannot parse VIZ_SERVER tag!" && tokens.size()==2);
      assert("pre: empty VIS_ZERGER tag" && !tokens[1].empty());

      this->Server = tokens[1];
      }
    else if(tokens[0] == "VIZ_PORT" )
      {
      assert("pre: cannot parse VIZ_PORT" && tokens.size()==2);
      assert("pre: empty VIZ_PORT tag" && !tokens[1].empty());

      this->Port = std::atoi(tokens[1].c_str());
      }
    else if(tokens[0] == "VIZ_FREQUENCY")
      {
      assert("pre: cannot parse VIZ_FREQUENCY" && tokens.size()==2);
      assert("pre: empty VIZ_FREQUENCY tag" && !tokens[1].empty());

      this->VisualizationFrequency = std::atoi(tokens[1].c_str());
      }
    else if( tokens[0] == "ANALYSISTOOL" )
      {
      assert("pre: cannot parse ANALYSISTOOL" && (tokens.size()==3));
      std::string toolName = tokens[1];
      bool toolStatus      = this->Parsebool(tokens[2]);
      this->RegisterAnalysisTool(toolName,toolStatus);
      }
    else if( tokens[0] == "#" && tokens.size()==3)
      {

      if(tokens[2] == "SECTION" )
        {
        currentAnalysisTool = tokens[1];
        if(!this->ToolExists(currentAnalysisTool))
          {
          std::cerr << "TOOL " << currentAnalysisTool << " is not registered!";
          std::cerr << std::endl;
          std::cerr << "FILE: " << __FILE__ << std::endl;
          std::cerr << "LINE: " << __LINE__ << std::endl;
          return false;
          } // END if ToolExists
        }// END if SECTION

      }
    else if( (currentAnalysisTool != "") && (tokens.size() > 1) )
      {
      if( !this->ToolExists(currentAnalysisTool) )
        {
        std::cerr << "ERROR: " << currentAnalysisTool << " tool doesn't exist!";
        std::cerr << "\nFILE: " << __FILE__ << std::endl;
        std::cerr << "LINE: " << __LINE__ << std::endl;
        return false;
        }

      if( tokens[0].substr(0,1) != "#" )
        {
        if( tokens.size() == 2 )
          {
          this->ToolToDictionary[currentAnalysisTool][tokens[0]]=tokens[1];
          }
        else
          {
          std::ostringstream oss;
          oss.clear();
          for( int i=1; i < tokens.size(); ++i )
            {
            oss << tokens[i] << " ";
            } // END for all value tokens
          this->ToolToDictionary[currentAnalysisTool][tokens[0]]=oss.str();
          }
        } // END if not a comment

      }

   } // END while

  return( detectedVersion );
}

//------------------------------------------------------------------------------
void CosmologyToolsConfiguration::WriteConfiguration(std::string file)
{
  std::ofstream ofs;
  ofs.open( file.c_str() );
  assert("pre: cannot open output file!" && ofs.is_open() );

  ofs << "#-------------------------------------------------------------------------------\n";
  ofs << "#  COSMOLOGY TOOLS CONFIGURATION\n";
  ofs << "#-------------------------------------------------------------------------------\n\n";

  ofs << "# Set the version of the configuration file (for backwards compatibility etc.)\n";
  ofs << "VERSION " << this->Version << std::endl;

  ofs << "# Visualization Parameters\n";
  ofs << "VISUALIZATION ";
  if( this->Visualization )
    {
    ofs << "YES\n";
    }
  else
    {
    ofs << "NO\n";
    }

  ofs << "SERVER " << this->Server << std::endl;
  ofs << "PORT "   << this->Port   << std::endl << std::endl;

  ofs << "## Frequency at which to update visualization\n";
  ofs << "VIZ_FREQUENCY " << this->VisualizationFrequency << std::endl;
  ofs << std::endl;

  std::map<std::string,bool>::iterator toolStatusIter;
  toolStatusIter = this->ToolStatus.begin();
  for(; toolStatusIter != this->ToolStatus.end(); ++toolStatusIter )
    {
    ofs << toolStatusIter->first << " ";
    if( toolStatusIter->second )
      {
      ofs << "YES\n";
      }
    else
      {
      ofs << "NO\n";
      }
    } // END for all tool statuses
  ofs << std::endl << std::endl;

  std::map<std::string,Dictionary>::iterator toolIterator;
  toolIterator = this->ToolToDictionary.begin();
  for( ;toolIterator != this->ToolToDictionary.end(); ++toolIterator)
    {
    ofs << "#-------------------------------------------------------------------------------\n";
    ofs << "# " << toolIterator->first << " SECTION" << std::endl;
    ofs << "#-------------------------------------------------------------------------------\n";

    Dictionary *toolParams = &toolIterator->second;
    assert("pre: tool parameters is NULL!" && (toolParams != NULL) );

    DictionaryIterator paramIter = toolParams->begin();
    for( ;paramIter != toolParams->end(); ++paramIter)
      {
      ofs << paramIter->first << " " << paramIter->second << std::endl;
      } // END for all parameter iterators

    ofs << std::endl;
    } // END for all tools

  ofs.close();
}

//------------------------------------------------------------------------------
bool CosmologyToolsConfiguration::ParseFile()
{
  assert("pre: No MPI Communicator was supplied!" &&
           (this->Communicator != MPI_COMM_NULL) );

  this->Clear();

  int rank = -1;
  MPI_Comm_rank(this->Communicator,&rank);

  int rawDataSize = 0;
  char *buffer    = NULL;
  switch( rank )
    {
    case 0:
      this->ReadInConfigurationFile();
      assert("pre: invalid configuration!" && !this->RawConfigData.empty());
      rawDataSize = static_cast<int>(this->RawConfigData.size()+1);
      MPI_Bcast(&rawDataSize,1,MPI_INT,0,this->Communicator);
      MPI_Bcast(const_cast<char*>(this->RawConfigData.c_str()),rawDataSize,
          MPI_CHAR,0,this->Communicator);
      break;
    default:
      MPI_Bcast(&rawDataSize,1,MPI_INT,0,this->Communicator);

      buffer = new char[rawDataSize];
      MPI_Bcast(buffer,rawDataSize,MPI_CHAR,0,this->Communicator);
      this->RawConfigData = std::string( buffer );
      delete [] buffer;
    } // END switch

  return( this->ParseRawData() );
}

//------------------------------------------------------------------------------
int CosmologyToolsConfiguration::GetNumberOfAnalysisTools()
{
  return(static_cast<int>(this->ToolNames.size()));
}

//------------------------------------------------------------------------------
std::string CosmologyToolsConfiguration::GetToolInstanceName(
    const int idx)
{
  std::string key = this->GetToolName(idx);
  return( this->GetToolClassInstance(key) );
}

//------------------------------------------------------------------------------
std::string CosmologyToolsConfiguration::GetToolClassInstance(
      std::string name)
{
  std::string instance ="";
  if(this->ToolToDictionary[name].find("INSTANCE_NAME") ==
      this->ToolToDictionary[name].end() )
    {
    instance = "NOT-FOUND";
    } // END if
  else
    {
    instance = this->ToolToDictionary[name]["INSTANCE_NAME"];
    }
  return( instance );
}

//------------------------------------------------------------------------------
std::string CosmologyToolsConfiguration::GetToolName(const int idx)
{
  assert("pre: index is out-of-bounds!" &&
          (idx >= 0) && (idx < this->GetNumberOfAnalysisTools()) );
  return( this->ToolNames[idx] );
}

//------------------------------------------------------------------------------
void CosmologyToolsConfiguration::GetToolParameters(
        std::string toolName, Dictionary &parameters)
{
  if(this->ToolToDictionary.find(toolName) != this->ToolToDictionary.end())
    {
    parameters = this->ToolToDictionary[toolName];
    }
}

//------------------------------------------------------------------------------
bool CosmologyToolsConfiguration::ToolExists(std::string toolName)
{
  bool status = false;
  if(this->ToolStatus.find(toolName) != this->ToolStatus.end() )
    {
    status = true;
    }
  else
    {
    status = false;
    }
  return( status );
}

//------------------------------------------------------------------------------
bool CosmologyToolsConfiguration::IsToolEnabled(std::string toolName)
{
  bool status = false;
  if( this->ToolExists( toolName ) )
    {
    status = this->ToolStatus[ toolName ];
    }
  else
    {
    status = false;
    }
  return( status );
}

} /* namespace cosmotk */
