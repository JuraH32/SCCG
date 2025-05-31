#ifndef POSITION_H
#define POSITION_H

#include <cstddef>
#include <limits>

struct Position {
    int startInReference;
    int endInReference;
    int startInTarget;
    int endInTarget;
    
    // Constructor allows for any of the fields to be initialized
    // with -1 to indicate "not set" or "invalid" state.
    explicit Position(
        int startInRef = -1,
        int endInRef = -1,
        int startInTar = -1,
        int endInTar = -1
    ) : 
        startInReference(startInRef),
        endInReference(endInRef),
        startInTarget(startInTar),
        endInTarget(endInTar) {}
    
    // Convenience getters and setters with more descriptive names
    [[nodiscard]] int getStartInReference() const { return startInReference; }
    void setStartInReference(int pos) { startInReference = pos; }

    [[nodiscard]] int getEndInReference() const { return endInReference; }
    void setEndInReference(int pos) { endInReference = pos; }

    [[nodiscard]] int getStartInTarget() const { return startInTarget; }
    void setStartInTarget(int pos) { startInTarget = pos; }

    [[nodiscard]] int getEndInTarget() const { return endInTarget; }
    void setEndInTarget(int pos) { endInTarget = pos; }
    
    // Check if a particular position is set
    [[nodiscard]] bool isStartInReferenceSet() const { return startInReference != -1; }
    [[nodiscard]] bool isEndInReferenceSet() const { return endInReference != -1; }
    [[nodiscard]] bool isStartInTargetSet() const { return startInTarget != -1; }
    [[nodiscard]] bool isEndInTargetSet() const { return endInTarget != -1; }

    [[nodiscard]] int length() const {
        if (isStartInTargetSet() && isEndInTargetSet()) {
            return endInTarget - startInTarget + 1;
        } else if (isStartInReferenceSet() && isEndInReferenceSet()) {
            return endInReference - startInReference + 1;
        } else {
            return -1;
        }

    }
};

#endif // POSITION_H
