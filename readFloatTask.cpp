#include "readFloatTask.h"

string readFloatTask::m_leftstr;
mutex readFloatTask::m_readmtx;
mutex readFloatTask::m_mergemtx;

readFloatTask::readFloatTask(/* args */)
{
    m_pbuf = new char[READLEN];
}

readFloatTask::~readFloatTask()
{
    delete[] m_pbuf;
}

void readFloatTask::doTask()
{
    bool brun = true;
    while (brun)
    {
        brun = readbuf();
        getInfo();
    }

    mergeinfo();
    printf("read task end\n");
}

bool readFloatTask::readbuf()
{
    lock_guard<mutex> lock(m_readmtx);
    char *pbuf = m_pbuf;
    int leftlen = m_leftstr.length();
    memcpy(pbuf, m_leftstr.c_str(), leftlen);
    m_leftstr.clear();
    pbuf += leftlen;
    int readlen = READLEN-leftlen;
    int reallen = gzread(cgef3dParam::GetInstance()->m_infile, pbuf, readlen);

    m_buflen = reallen;
    if(reallen == readlen)
    {
        cuttail(m_pbuf);
    }
    else
    {
        if (m_buflen != 0)
            m_buflen += leftlen;
        return false;
    }
    return true;
}

int readFloatTask::cuttail(char *pbuf)
{
    int i = READLEN-1;
    for(;i>0;i--)
    {
        if(pbuf[i] == '\n')
        {
            break;
        }
    }

    m_buflen = i+1;
    m_leftstr.append(&pbuf[m_buflen], READLEN-m_buflen);
    return 0;
}

int readFloatTask::getInfo()
{
    uint32_t i = 0, k = 0, celllabel=0, len = 0;
    char *ptr = m_pbuf;

    string gname;
    float x=0.0,y=0.0,umi=0.0;
    for(;i<m_buflen;i++)
    {
        if(m_pbuf[i] == '\t' || m_pbuf[i] == '\n')
        {
            switch (k)
            {
            case 0:
                len = &m_pbuf[i]-ptr;
                gname.clear();
                gname.append(ptr, len);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 1:
                x = atof(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 2:
                y = atof(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 3:
                umi = atof(ptr);
                k++;
                ptr = &m_pbuf[i+1];
                break;
            case 4:
                k = 0;
                celllabel = atoi(ptr);
                ptr = &m_pbuf[i+1];

                if(celllabel > 0)
                {
                    if(m_map_cell.find(celllabel) == m_map_cell.end())
                    {
                        cgef3d_cell *ptr = new cgef3d_cell();
                        m_map_cell.emplace(celllabel, ptr);
                    }
                    m_map_cell[celllabel]->add(x,y,umi);

                    if(m_map_gene.find(gname) == m_map_gene.end())
                    {
                        cgef3d_gene *gptr = new cgef3d_gene();
                        m_map_gene.emplace(gname, gptr);
                    }
                    m_map_gene[gname]->add(celllabel, umi);
                }
                break;
            default:
                break;
            }
        }
    }

    return m_map_cell.size();
}

int readFloatTask::mergeinfo()
{
    lock_guard<mutex> lock(m_mergemtx);

    auto itor = m_map_cell.begin();
    auto &tmap_cell = cgef3dParam::GetInstance()->m_map_cell;
    for(;itor!=m_map_cell.end();itor++)
    {
        if(tmap_cell.find(itor->first) != tmap_cell.end())
        {
            tmap_cell[itor->first]->merge(*(itor->second));
            delete itor->second;
        }
        else
        {
            tmap_cell.emplace(itor->first, itor->second);
        }
    }

    auto itor_g = m_map_gene.begin();
    auto &tmp_gene = cgef3dParam::GetInstance()->m_map_gene;
    for(;itor_g!=m_map_gene.end();itor_g++)
    {
        if(tmp_gene.find(itor_g->first) != tmp_gene.end())
        {
            tmp_gene[itor_g->first]->merge(*(itor_g->second));
            delete itor_g->second;
        }
        else
        {
            tmp_gene.emplace(itor_g->first, itor_g->second);
        }
    }

    // auto itor_c = m_map_ctype.begin();
    // auto &tmp_ctype = cgef3dParam::GetInstance()->m_map_ctype;
    // for(;itor_c!=m_map_ctype.end();itor_c++)
    // {
    //     if(tmp_ctype.find(itor_c->first) != tmp_ctype.end())
    //     {
    //         tmp_ctype[itor_g->first]->merge(*(itor_c->second));
    //         delete itor_c->second;
    //     }
    //     else
    //     {
    //         tmp_ctype.emplace(itor_c->first, itor_c->second);
    //     }
    // }

    return 0;
}

// int readFloatTask_5::getInfo()
// {
//     uint32_t i = 0, k = 0, celllabel=0, len = 0;
//     char *ptr = m_pbuf;

//     string gname;
//     float x=0.0,y=0.0,umi=0.0;
//     for(;i<m_buflen;i++)
//     {
//         if(m_pbuf[i] == '\t' || m_pbuf[i] == '\n')
//         {
//             switch (k)
//             {
//             case 0:
//                 len = &m_pbuf[i]-ptr;
//                 gname.clear();
//                 gname.append(ptr, len);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 1:
//                 x = atof(ptr);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 2:
//                 y = atof(ptr);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 3:
//                 umi = atof(ptr);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 4:
//                 k = 0;
//                 celllabel = atoi(ptr);
//                 ptr = &m_pbuf[i+1];

//                 if(celllabel > 0)
//                 {
//                     if(m_map_cell.find(celllabel) == m_map_cell.end())
//                     {
//                         cgef3d_cell *ptr = new cgef3d_cell();
//                         m_map_cell.emplace(celllabel, ptr);
//                     }
//                     m_map_cell[celllabel]->add(x,y,umi);

//                     if(m_map_gene.find(gname) == m_map_gene.end())
//                     {
//                         cgef3d_gene *gptr = new cgef3d_gene();
//                         m_map_gene.emplace(gname, gptr);
//                     }
//                     m_map_gene[gname]->add(celllabel, umi);
//                 }
//                 break;
//             default:
//                 break;
//             }
//         }
//     }

//     return m_map_cell.size();
// }

// int readFloatTask_6::getInfo()
// {
//     uint32_t i = 0, k = 0, celllabel=0, len = 0;
//     char *ptr = m_pbuf;

//     string gname, ctype;
//     float x=0.0,y=0.0,umi=0.0;
//     for(;i<m_buflen;i++)
//     {
//         if(m_pbuf[i] == '\t' || m_pbuf[i] == '\n')
//         {
//             switch (k)
//             {
//             case 0:
//                 len = &m_pbuf[i]-ptr;
//                 gname.clear();
//                 gname.append(ptr, len);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 1:
//                 x = atof(ptr);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 2:
//                 y = atof(ptr);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 3:
//                 umi = atof(ptr);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 4:
//                 celllabel = atoi(ptr);
//                 k++;
//                 ptr = &m_pbuf[i+1];
//                 break;
//             case 5:
//                 k = 0;
//                 if(celllabel > 0)
//                 {
//                     if(m_map_cell.find(celllabel) == m_map_cell.end())
//                     {
//                         cgef3d_cell *ptr = new cgef3d_cell();
//                         m_map_cell.emplace(celllabel, ptr);
//                     }
//                     m_map_cell[celllabel]->add(x,y,umi);

//                     if(m_map_gene.find(gname) == m_map_gene.end())
//                     {
//                         cgef3d_gene *gptr = new cgef3d_gene();
//                         m_map_gene.emplace(gname, gptr);
//                     }
//                     m_map_gene[gname]->add(celllabel, umi);

//                     ctype.clear();
//                     ctype.append( ptr, &m_pbuf[i]-ptr);
//                     if(m_map_ctype.find(ctype) == m_map_ctype.end())
//                     {
//                         cgef3d_ctype *cptr = new cgef3d_ctype();
//                         m_map_ctype.emplace(ctype, cptr);
//                     }
//                     m_map_ctype[ctype]->add(celllabel);
//                 }
//                 ptr = &m_pbuf[i+1];
//                 break;
//             default:
//                 break;
//             }
//         }
//     }

//     return m_map_cell.size();
// }