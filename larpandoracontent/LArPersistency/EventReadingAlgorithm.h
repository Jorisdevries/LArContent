/**
 *  @file   larpandoracontent/LArPersistency/EventReadingAlgorithm.h
 * 
 *  @brief  Header file for the event reading algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_EVENT_READING_ALGORITHM_H
#define LAR_EVENT_READING_ALGORITHM_H 1

#include "Pandora/ExternallyConfiguredAlgorithm.h"

#include "Pandora/PandoraInputTypes.h"

#include "Persistency/PandoraIO.h"

namespace pandora {class FileReader;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  EventReadingAlgorithm class
 */
class EventReadingAlgorithm : public pandora::ExternallyConfiguredAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventReadingAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventReadingAlgorithm();

    /**
     *  @brief  External event reading parameters class
     */
    class ExternalEventReadingParameters : public ExternalParameters
    {
    public:
        std::string             m_geometryFileName;             ///< Name of the file containing geometry information
        std::string             m_eventFileName;                ///< Name of the file containing event information
        pandora::InputUInt      m_skipToEvent;                  ///< Index of first event to consider in input file
    };

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode Run();

    /**
     *  @brief  Analyze a provided file name to extract the file type/extension
     *
     *  @param  fileName the file name
     *
     *  @return the file type
     */
    pandora::FileType GetFileType(const std::string &fileName) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::FileType           m_geometryFileType;             ///< The geometry file type
    pandora::FileType           m_eventFileType;                ///< The event file type

    std::string                 m_geometryFileName;             ///< Name of the file containing geometry information
    std::string                 m_eventFileName;                ///< Name of the file containing event information
    unsigned int                m_skipToEvent;                  ///< Index of first event to consider in input file

    bool                        m_useLArCaloHits;               ///< Whether to read lar calo hits, or standard pandora calo hits
    bool                        m_useLArMCParticles;            ///< Whether to read lar mc particles, or standard pandora mc particles

    pandora::FileReader        *m_pEventFileReader;             ///< Address of the event file reader
};

} // namespace lar_content

#endif // #ifndef LAR_EVENT_READING_ALGORITHM_H
