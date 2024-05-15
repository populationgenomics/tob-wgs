from collections import defaultdict

from metamist.apis import AssayApi
from metamist.graphql import gql, query
from metamist.models import AssayUpsert

QUERY = gql(
    """
        query MyQuery {
            project(name: "tob-wgs-dev") {
                participants {
                id
                samples {
                    assays {
                        id
                    }
                }
                }
            }
        }
    """,
)


def get_data_metamist():
    participants = query(QUERY)['project']['participants']
    list_assays = defaultdict(list)

    for participant in participants:
        for sample in participant['samples']:
            list_assays[participant['id']].append(sample['assays'][0]['id'])

    return list_assays


def upsert_tube_ids(list_participants: defaultdict):
    print(
        'list_participants is hard coded in the function. Ensure dummy_tube keys are representative of participants in list',
    )
    print(list_participants)
    aapi = AssayApi()
    # Create dictionary containing dummy tube IDs
    dummy_tubes = {
        1368: 'FDS1000000000',
        1369: 'FDS1000000001',
        1370: 'FDS1000000002',
        1371: 'FDS1000000003',
        1372: 'FDS1000000004',
    }

    for key, tube in dummy_tubes.items():
        new_assay = AssayUpsert()
        new_assay['id'] = key
        new_assay['meta'] = {}
        new_assay['meta']['kgcc_tube_id'] = tube
        aapi.update_assay(assay_upsert=new_assay)


if __name__ == ('__main__'):
    list_participants = get_data_metamist()
    upsert_tube_ids(list_participants)
